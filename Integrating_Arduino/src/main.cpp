#include <Arduino.h>
#include "imu.h"
#include "madgwick.h"
#include "orientation.h"
#include "telemetry.h"
#include "motion_pipeline.h"

// =====================================================
//                      SETTINGS
// =====================================================

// --------- Rates ---------
static const float SERIAL_HZ = 50.0f;    // CSV rate
static const float BLE_HZ    = 50.0f;    // BLE rate
static const float LOOP_DT_MAX = 0.05f;  // clamp dt to 50 ms
static const float LOOP_DT_MIN = 1e-4f;

// --------- Madgwick ---------
static const float BETA_NOM     = 0.10f;
static const float BETA_SETTLE  = 0.30f;
static const float ROT_TRIG_RAD = 0.35f;
static const float SETTLE_S     = 0.60f;

// --------- Boot gyro bias calib ---------
static const uint32_t GYRO_BIAS_CAL_MS = 2500;
static bool  g_gyroBiasReady = false;
static float g_gxBias=0.0f, g_gyBias=0.0f, g_gzBias=0.0f;

// --------- Runtime toggles ---------
static bool g_enableSerialCSV = true;
static bool g_enableBetaBoost = true;
static bool g_enableZUPT      = false;

// --------- Presets ---------
enum FilterPreset : uint8_t { PRESET_NORMAL=0, PRESET_AGGRESSIVE=1, PRESET_LIGHT=2 };
static FilterPreset g_preset = PRESET_NORMAL;

// =====================================================
//                  UTILS / HELPERS
// =====================================================


static void printHelp() {
  Serial.println("Commands:");
  Serial.println("  r  -> reset pipeline (pos/vel=0, timers)");
  Serial.println("  z  -> toggle ZUPT");
  Serial.println("  b  -> toggle beta-boost");
  Serial.println("  p  -> toggle Serial CSV output");
  Serial.println("  t  -> print status line");
  Serial.println("  0  -> preset NORMAL");
  Serial.println("  1  -> preset AGGRESSIVE");
  Serial.println("  2  -> preset LIGHT");
}

static void calibrateGyroBias() {
  uint32_t t0 = millis();
  uint32_t n = 0;
  double sx=0, sy=0, sz=0;

  Serial.println("Calibrating gyro bias... keep still");
  while (millis() - t0 < GYRO_BIAS_CAL_MS) {
    float ax,ay,az,gx,gy,gz,mx,my,mz;
    imuRead(ax,ay,az,gx,gy,gz,mx,my,mz);
    sx += gx; sy += gy; sz += gz;
    n++;
    delay(10);
  }

  if (n > 0) {
    g_gxBias = (float)(sx / (double)n);
    g_gyBias = (float)(sy / (double)n);
    g_gzBias = (float)(sz / (double)n);
    g_gyroBiasReady = true;
  }

  Serial.print("Gyro bias (rad/s): ");
  Serial.print(g_gxBias, 6); Serial.print(", ");
  Serial.print(g_gyBias, 6); Serial.print(", ");
  Serial.println(g_gzBias, 6);
}

static PipelineConfig makePreset(FilterPreset p) {
  PipelineConfig cfg;

  // defaults (NORMAL)
  cfg.acc_lpf_hz = 12.0f;
  cfg.gyro_lpf_hz = 12.0f;
  cfg.hampel_nsigma = 4.0f;
  cfg.acc_bias_tau_s = 20.0f;

  cfg.still.gyro_still_rad_s = 0.08f;
  cfg.still.amag_still_err_g = 0.05f;
  cfg.still.hold_s = 0.35f;
  cfg.still.enable_zupt = true;

  cfg.integ.dt_min_s = 1e-4f;
  cfg.integ.dt_max_s = 0.20f;
  cfg.integ.acc_hpf_hz = 0.10f;
  cfg.integ.vel_leak_tau_s = 6.0f;
  cfg.integ.jerk_max_mps3 = 80.0f;

  if (p == PRESET_AGGRESSIVE) {
    cfg.acc_lpf_hz = 8.0f;
    cfg.gyro_lpf_hz = 8.0f;
    cfg.hampel_nsigma = 3.5f;
    cfg.acc_bias_tau_s = 12.0f;

    cfg.still.gyro_still_rad_s = 0.06f;
    cfg.still.amag_still_err_g = 0.04f;
    cfg.still.hold_s = 0.25f;

    cfg.integ.acc_hpf_hz = 0.15f;
    cfg.integ.vel_leak_tau_s = 4.0f;
    cfg.integ.jerk_max_mps3 = 60.0f;
  } else if (p == PRESET_LIGHT) {
    cfg.acc_lpf_hz = 18.0f;
    cfg.gyro_lpf_hz = 18.0f;
    cfg.hampel_nsigma = 5.0f;
    cfg.acc_bias_tau_s = 30.0f;

    cfg.still.gyro_still_rad_s = 0.10f;
    cfg.still.amag_still_err_g = 0.06f;
    cfg.still.hold_s = 0.45f;

    cfg.integ.acc_hpf_hz = 0.07f;
    cfg.integ.vel_leak_tau_s = 8.0f;
    cfg.integ.jerk_max_mps3 = 90.0f;
  }

  return cfg;
}

// =====================================================
//                        MAIN
// =====================================================

int main(void) {
  init();

  Serial.begin(115200);
  delay(1200);

  Serial.println("=== Glider IMU + Motion Pipeline ===");

  if (!imuBegin()) {
    Serial.println("IMU init failed.");
    while (1) delay(1000);
  }

  // BLE init (2 characteristics: IMU14 + NAV10)
  bool bleOk = telemetryBegin();
  Serial.print("BLE begin: ");
  Serial.println(bleOk ? "OK" : "FAIL");

  printHelp();
  calibrateGyroBias();

  // --------- Serial CSV header ---------
  Serial.println(
    "t_s,"
    "qw,qx,qy,qz,"
    "axg,ayg,azg,"
    "gx,gy,gz,"
    "axlin_g,aylin_g,azlin_g,"
    "axw,ayw,azw,"
    "vx,vy,vz,"
    "px,py,pz,"
    "roll,pitch,yaw,"
    "still,zupt,reject,jerk,"
    "dt,dt_clamped,ble"
  );

  // --------- Madgwick ---------
  MadgwickAHRS ahrs(50.0f, BETA_NOM);
  float settle_timer_s = 0.0f;

  // --------- Pipeline ---------
  MotionPipeline pipe;
  PipelineConfig pcfg = makePreset(g_preset);
  pcfg.still.enable_zupt = g_enableZUPT;
  pipe.setConfig(pcfg);
  pipe.reset();

  // --------- Timers / schedulers ---------
  uint32_t t0_us = micros();
  uint32_t last_us = t0_us;

  const float serial_period_s = 1.0f / SERIAL_HZ;
  const float ble_period_s    = 1.0f / BLE_HZ;

  float serial_acc_s = 0.0f;
  float ble_acc_s    = 0.0f;

  // --------- Debug counters ---------
  uint32_t dt_glitch_count = 0;
  uint32_t rejected_count  = 0;
  uint32_t zupt_count      = 0;

  // --------- Last out (for status) ---------
  PipelineOut out{};
  float last_dt = 0.0f;
  float last_dt_clamped = 0.0f;

  while (1) {
    // =================================================
    // 1) Read Serial commands
    // =================================================
    while (Serial.available()) {
      char c = (char)Serial.read();
      if (c == '\n' || c == '\r') continue;

      if (c == 'h' || c == '?') {
        printHelp();
      } else if (c == 'r' || c == 'R') {
        pipe.reset();
        settle_timer_s = 0.0f;
        serial_acc_s = 0.0f;
        ble_acc_s = 0.0f;
        dt_glitch_count = 0;
        rejected_count  = 0;
        zupt_count      = 0;
        Serial.println("RESET OK");
      } else if (c == 'z' || c == 'Z') {
        g_enableZUPT = !g_enableZUPT;
        pcfg.still.enable_zupt = g_enableZUPT;
        pipe.setConfig(pcfg);
        Serial.print("ZUPT: ");
        Serial.println(g_enableZUPT ? "ON" : "OFF");
      } else if (c == 'b' || c == 'B') {
        g_enableBetaBoost = !g_enableBetaBoost;
        Serial.print("Beta boost: ");
        Serial.println(g_enableBetaBoost ? "ON" : "OFF");
      } else if (c == 'p' || c == 'P') {
        g_enableSerialCSV = !g_enableSerialCSV;
        Serial.print("Serial CSV: ");
        Serial.println(g_enableSerialCSV ? "ON" : "OFF");
      } else if (c == 't' || c == 'T') {
        Serial.print("preset="); Serial.print((int)g_preset);
        Serial.print("  ZUPT="); Serial.print(g_enableZUPT ? "1" : "0");
        Serial.print("  betaBoost="); Serial.print(g_enableBetaBoost ? "1" : "0");
        Serial.print("  BLE="); Serial.print(telemetryIsConnected() ? "1" : "0");
        Serial.print("  rejects="); Serial.print(rejected_count);
        Serial.print("  zupt="); Serial.print(zupt_count);
        Serial.print("  dtGlitch="); Serial.println(dt_glitch_count);
      } else if (c == '0' || c == '1' || c == '2') {
        g_preset = (c == '0') ? PRESET_NORMAL : (c == '1') ? PRESET_AGGRESSIVE : PRESET_LIGHT;
        pcfg = makePreset(g_preset);
        pcfg.still.enable_zupt = g_enableZUPT;
        pipe.setConfig(pcfg);
        Serial.print("Preset set to ");
        Serial.println((int)g_preset);
      }
    }

    telemetryPoll();
    // =================================================
    // 2) Timing (dt)
    // =================================================
    uint32_t now_us = micros();
    float dt = (now_us - last_us) * 1e-6f;
    last_us = now_us;

    last_dt = dt;

    bool dt_clamped_flag = false;
    if (!isfinite(dt) || dt <= 0.0f) {
      dt = LOOP_DT_MIN;
      dt_clamped_flag = true;
      dt_glitch_count++;
    }
    if (dt > LOOP_DT_MAX) {
      dt = LOOP_DT_MAX;
      dt_clamped_flag = true;
      dt_glitch_count++;
    }
    if (dt < LOOP_DT_MIN) {
      dt = LOOP_DT_MIN;
      dt_clamped_flag = true;
      dt_glitch_count++;
    }

    last_dt_clamped = dt;

    serial_acc_s += dt;
    ble_acc_s    += dt;

    // =================================================
    // 3) Read IMU
    // =================================================
    float ax, ay, az, gx, gy, gz, mx, my, mz;
    imuRead(ax,ay,az,gx,gy,gz,mx,my,mz);

    if (g_gyroBiasReady) { gx -= g_gxBias; gy -= g_gyBias; gz -= g_gzBias; }

    // =================================================
    // 4) Madgwick beta boost
    // =================================================
    float beta = BETA_NOM;
    if (g_enableBetaBoost) {
      float wmag = sqrtf(gx*gx + gy*gy + gz*gz);
      if (wmag > ROT_TRIG_RAD) settle_timer_s = SETTLE_S;
      settle_timer_s = max(0.0f, settle_timer_s - dt);
      beta = (settle_timer_s > 0.0f) ? BETA_SETTLE : BETA_NOM;
    }
    ahrs.beta = beta;

    // Orientation update
    ahrs.updateIMU(gx, gy, gz, ax, ay, az);
    float qw = ahrs.q0, qx = ahrs.q1, qy = ahrs.q2, qz = ahrs.q3;

    // =================================================
    // 5) Motion pipeline
    // =================================================
    pipe.update(dt, qw,qx,qy,qz, {ax,ay,az}, {gx,gy,gz}, out);

    if (out.rejected) rejected_count++;
    if (out.zupt_applied) zupt_count++;

    float t_s = (now_us - t0_us) * 1e-6f;

    // =================================================
    // 6) Serial CSV (throttled)
    // =================================================
    if (g_enableSerialCSV && serial_acc_s >= serial_period_s) {
      serial_acc_s = 0.0f;

      Serial.print(t_s, 6); Serial.print(',');

      Serial.print(out.qw,6); Serial.print(',');
      Serial.print(out.qx,6); Serial.print(',');
      Serial.print(out.qy,6); Serial.print(',');
      Serial.print(out.qz,6); Serial.print(',');

      Serial.print(out.a_g_body.x,6); Serial.print(',');
      Serial.print(out.a_g_body.y,6); Serial.print(',');
      Serial.print(out.a_g_body.z,6); Serial.print(',');

      Serial.print(out.w_rad_body.x,6); Serial.print(',');
      Serial.print(out.w_rad_body.y,6); Serial.print(',');
      Serial.print(out.w_rad_body.z,6); Serial.print(',');

      Serial.print(out.a_lin_g_body.x,6); Serial.print(',');
      Serial.print(out.a_lin_g_body.y,6); Serial.print(',');
      Serial.print(out.a_lin_g_body.z,6); Serial.print(',');

      Serial.print(out.a_lin_mps2_world.x,6); Serial.print(',');
      Serial.print(out.a_lin_mps2_world.y,6); Serial.print(',');
      Serial.print(out.a_lin_mps2_world.z,6); Serial.print(',');

      Serial.print(out.v_w.x,6); Serial.print(',');
      Serial.print(out.v_w.y,6); Serial.print(',');
      Serial.print(out.v_w.z,6); Serial.print(',');

      Serial.print(out.p_w.x,6); Serial.print(',');
      Serial.print(out.p_w.y,6); Serial.print(',');
      Serial.print(out.p_w.z,6); Serial.print(',');

      Serial.print(out.roll_deg,3); Serial.print(',');
      Serial.print(out.pitch_deg,3); Serial.print(',');
      Serial.print(out.yaw_deg,3); Serial.print(',');

      Serial.print(out.still ? 1 : 0); Serial.print(',');
      Serial.print(out.zupt_applied ? 1 : 0); Serial.print(',');
      Serial.print(out.rejected ? 1 : 0); Serial.print(',');
      Serial.print(out.jerk,3); Serial.print(',');

      Serial.print(last_dt,6); Serial.print(',');
      Serial.print(dt_clamped_flag ? 1 : 0); Serial.print(',');
      Serial.print(telemetryIsConnected() ? 1 : 0);

      Serial.println();
    }

    // =================================================
    // 7) BLE (throttled)
    // =================================================
    if (telemetryIsConnected() && ble_acc_s >= ble_period_s) {
      ble_acc_s = 0.0f;

      // IMU14: 14 floats (56B)
      float imu14[14] = {
        t_s,
        out.qw, out.qx, out.qy, out.qz,
        out.a_lin_mps2_world.x, out.a_lin_mps2_world.y, out.a_lin_mps2_world.z,
        out.w_rad_body.x, out.w_rad_body.y, out.w_rad_body.z,
        out.roll_deg, out.pitch_deg, out.yaw_deg
      };
      telemetryUpdateIMU14(imu14);

      // NAV10: 10 floats (40B)
      float nav10[10] = {
        out.v_w.x, out.v_w.y, out.v_w.z,
        out.p_w.x, out.p_w.y, out.p_w.z,
        out.still ? 1.0f : 0.0f,
        out.zupt_applied ? 1.0f : 0.0f,
        out.rejected ? 1.0f : 0.0f,
        out.jerk
      };
      telemetryUpdateNAV10(nav10);
    }

    // Small yield
    delay(1);
  }

  return 0;
}
