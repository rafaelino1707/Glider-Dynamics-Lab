#include <Arduino.h>
#include <ArduinoBLE.h>
#include "imu.h"
#include "madgwick.h"
#include "orientation.h"
#include "telemetry.h"

// -----------------------------------------------------
// Operation Mode
// -----------------------------------------------------
#define ECO_MODE 2

// -----------------------------------------------------
// Madgwick tuning
// -----------------------------------------------------
static const float BETA_NOM        = 0.10f;
static const float BETA_DYNAMIC    = 0.06f;
static const float BETA_SETTLE     = 0.30f;
static const float BETA_SOFT_SCALE = 0.20f;

// amag thresholds (amag in g)
static const float AMAG_THR_SOFT = 0.15f;
static const float AMAG_THR_HARD = 0.25f;

// -----------------------------------------------------
// "Almost still" detector (for fast gravity convergence)
// -----------------------------------------------------
static const float AMAG_STILL_ERR_G   = 0.05f;   // |amag-1| < 0.05g
static const float GYRO_STILL_RAD_S   = 0.08f;   // |w| < 0.08 rad/s
static const float STILL_HOLD_S       = 0.35f;   // must hold for this time

// After a rotation ends, boost beta for some time to re-align gravity
static const float ROT_RATE_TRIG_RAD_S   = 0.35f;
static const float SETTLE_AFTER_ROT_S    = 0.60f;

// -----------------------------------------------------
// Simple sensor LPFs (pre-Madgwick)
// -----------------------------------------------------
static const float ACC_LPF_HZ  = 12.0f;
static const float GYRO_LPF_HZ = 12.0f;

// -----------------------------------------------------
// Gyro bias calibration (boot)
// -----------------------------------------------------
static const uint32_t GYRO_BIAS_CAL_MS = 2500;
static bool  gyro_bias_ready = false;
static float gx_bias = 0.0f, gy_bias = 0.0f, gz_bias = 0.0f;

// -----------------------------------------------------
// Target cadence (used as initial only; runtime becomes dynamic)
// -----------------------------------------------------
static const float SAMPLE_FREQ_TARGET = 30.0f;  // initial hint only

// -----------------------------------------------------
// Dynamic sampling (adaptive cadence)
// IMPORTANT: do NOT use still_flag to drop freq in flight.
// -----------------------------------------------------
static const float FREQ_MIN_HZ   = 25.0f;  // raise minimum so it never becomes "too slow"
static const float FREQ_MAX_HZ   = 60.0f;
static const float FREQ_SMOOTH_A = 0.18f;  // 0..1 (higher reacts faster)

// normalize "motion" scaling
static const float GYRO_NORM_FULL = 1.2f;  // rad/s ~ strong motion
static const float AMAG_ERR_FULL  = 0.30f; // g     ~ strong motion

static float    g_targetFreqHz   = SAMPLE_FREQ_TARGET;
static uint32_t g_samplePeriodUs = (uint32_t)(1000000.0f / SAMPLE_FREQ_TARGET);

// motion memory from previous cycle (0..1)
static float g_prevMotionMetric = 1.0f;

// -----------------------------------------------------
// Filter instance
// -----------------------------------------------------
MadgwickAHRS filter(SAMPLE_FREQ_TARGET, BETA_NOM);

// -----------------------------------------------------
// Logging window (ECO_MODE 1/2)
// -----------------------------------------------------
const uint32_t START_LOG_MS    = 120000UL;
const uint32_t LOG_DURATION_MS = 30000UL;

// -----------------------------------------------------
// RAM Logging
// -----------------------------------------------------
struct Sample {
  uint32_t t_ms;
  float qw, qx, qy, qz;
  float ax, ay, az;      // accel [g]
  float gx, gy, gz;      // gyro  [rad/s]
  float mx, my, mz;      // mag
};

constexpr size_t LOG_MAX_SAMPLES = 2000;
Sample g_logBuf[LOG_MAX_SAMPLES];
size_t g_logCount = 0;

static inline float clamp01(float x) { return (x < 0.0f) ? 0.0f : (x > 1.0f ? 1.0f : x); }
static inline float clampf(float x, float lo, float hi) { return (x < lo) ? lo : (x > hi ? hi : x); }

static inline void updateDynamicCadenceFromPrev() {
  // Use ONLY motion metric to set cadence (not still_flag).
  float desired = FREQ_MIN_HZ + (FREQ_MAX_HZ - FREQ_MIN_HZ) * g_prevMotionMetric;

  g_targetFreqHz = (1.0f - FREQ_SMOOTH_A) * g_targetFreqHz + (FREQ_SMOOTH_A) * desired;
  g_targetFreqHz = clampf(g_targetFreqHz, FREQ_MIN_HZ, FREQ_MAX_HZ);

  g_samplePeriodUs = (uint32_t)(1000000.0f / g_targetFreqHz);
  if (g_samplePeriodUs < 2000)   g_samplePeriodUs = 2000;    // 500 Hz cap safety
  if (g_samplePeriodUs > 200000) g_samplePeriodUs = 200000;  // 5 Hz cap safety
}

void dumpRamLogToSerial() {
  Serial.println("t_ms,qw,qx,qy,qz,ax,ay,az,gx,gy,gz,mx,my,mz");
  for (size_t i = 0; i < g_logCount; ++i) {
    const Sample &s = g_logBuf[i];
    Serial.print(s.t_ms); Serial.print(',');
    Serial.print(s.qw, 6); Serial.print(',');
    Serial.print(s.qx, 6); Serial.print(',');
    Serial.print(s.qy, 6); Serial.print(',');
    Serial.print(s.qz, 6); Serial.print(',');
    Serial.print(s.ax, 6); Serial.print(',');
    Serial.print(s.ay, 6); Serial.print(',');
    Serial.print(s.az, 6); Serial.print(',');
    Serial.print(s.gx, 6); Serial.print(',');
    Serial.print(s.gy, 6); Serial.print(',');
    Serial.print(s.gz, 6); Serial.print(',');
    Serial.print(s.mx, 6); Serial.print(',');
    Serial.print(s.my, 6); Serial.print(',');
    Serial.print(s.mz, 6);
    Serial.println();
  }
}

static void calibrateGyroBias() {
  uint32_t t0 = millis();
  uint32_t n = 0;
  double sx = 0.0, sy = 0.0, sz = 0.0;

  while (millis() - t0 < GYRO_BIAS_CAL_MS) {
    float ax, ay, az, gx, gy, gz, mx, my, mz;
    imuRead(ax, ay, az, gx, gy, gz, mx, my, mz);
    sx += gx; sy += gy; sz += gz;
    n++;
    delay(1000 / (int)SAMPLE_FREQ_TARGET);
  }

  if (n > 0) {
    gx_bias = (float)(sx / (double)n);
    gy_bias = (float)(sy / (double)n);
    gz_bias = (float)(sz / (double)n);
    gyro_bias_ready = true;
  }
}

// 1st-order LPF step
static inline float lpf1(float y, float x, float dt, float fc_hz) {
  if (fc_hz <= 0.0f || dt <= 0.0f) return x;
  float a = (2.0f * PI * fc_hz * dt) / (1.0f + 2.0f * PI * fc_hz * dt);
  return y + a * (x - y);
}

// Rotate vector v_body -> v_world using quaternion (qw,qx,qy,qz)
static inline void quatRotateBodyToWorld(float qw, float qx, float qy, float qz,
                                        float vx, float vy, float vz,
                                        float &ox, float &oy, float &oz)
{
  float qvx =  qy*vz - qz*vy;
  float qvy =  qz*vx - qx*vz;
  float qvz =  qx*vy - qy*vx;
  float qvw = -qx*vx - qy*vy - qz*vz;

  ox = vx + 2.0f*(qw*qvx + qvw*(-qx) + qvy*(-qz) - qvz*(-qy));
  oy = vy + 2.0f*(qw*qvy + qvw*(-qy) + qvz*(-qx) - qvx*(-qz));
  oz = vz + 2.0f*(qw*qvz + qvw*(-qz) + qvx*(-qy) - qvy*(-qx));
}

// Gravity vector in BODY frame from quaternion
static inline void gravityBodyFromQuat(float qw, float qx, float qy, float qz,
                                      float &gx_b, float &gy_b, float &gz_b)
{
  gx_b = 2.0f*(qx*qz - qw*qy);
  gy_b = 2.0f*(qw*qx + qy*qz);
  gz_b = qw*qw - qx*qx - qy*qy + qz*qz;
}

void setup() {
  Serial.begin(115200);
  delay(2000);

  Serial.println("Glider IMU + Madgwick (dt real + gyro bias + pre-LPF + settle beta + dynamic Hz)");

  if (!imuBegin()) {
    Serial.println("IMU init failed");
    while (1) delay(1000);
  }
  Serial.println("IMU OK");

  if (!telemetryBegin()) {
    Serial.println("BLE init failed");
    while (1) delay(1000);
  }
  Serial.println("BLE advertising");

  Serial.println("Calibrating gyro bias... keep still");
  calibrateGyroBias();
  if (gyro_bias_ready) {
    Serial.print("Gyro bias (rad/s): ");
    Serial.print(gx_bias, 6); Serial.print(", ");
    Serial.print(gy_bias, 6); Serial.print(", ");
    Serial.println(gz_bias, 6);
  }

  // CSV header (clean)
  Serial.println(
    "t_ms,"
    "qw,qx,qy,qz,"
    "ax_lin_w_mps2,ay_lin_w_mps2,az_lin_w_mps2,"
    "gx_rads,gy_rads,gz_rads,"
    "roll_deg,pitch_deg,yaw_deg,"
    "amag_g,gyro_norm,amag_err_g,"
    "beta_used,still_flag,settle_flag,"
    "targetFreq_Hz,sampleFreq_Hz"
  );
}

void loop() {
  if (Serial.available() > 0) {
    int c = Serial.read();
    if (c == 'D' || c == 'd') {
      Serial.println("=== RAM LOG DUMP BEGIN ===");
      dumpRamLogToSerial();
      Serial.println("=== RAM LOG DUMP END ===");
    }
  }

  BLE.poll();

  static uint32_t last_us = micros();
  uint32_t now_us = micros();
  uint32_t dt_us = now_us - last_us;

  // dynamic cadence gate (from previous motion metric)
  updateDynamicCadenceFromPrev();
  if (dt_us < g_samplePeriodUs) return;

  last_us = now_us;

  // real dt clamp
  if (dt_us < 500)    dt_us = 500;
  if (dt_us > 200000) dt_us = 200000;
  float dt_s = (float)dt_us * 1e-6f;

  float sampleFreq = 1e6f / (float)dt_us;
  filter.sampleFreq = sampleFreq;

  // read IMU
  float ax_g, ay_g, az_g;   // g
  float gx, gy, gz;         // rad/s
  float mx, my, mz;
  imuRead(ax_g, ay_g, az_g, gx, gy, gz, mx, my, mz);

  // gyro bias remove
  if (gyro_bias_ready) {
    gx -= gx_bias; gy -= gy_bias; gz -= gz_bias;
  }

  // pre-LPF
  static bool lpf_init = false;
  static float ax_f=0, ay_f=0, az_f=0;
  static float gx_f=0, gy_f=0, gz_f=0;

  if (!lpf_init) {
    ax_f=ax_g; ay_f=ay_g; az_f=az_g;
    gx_f=gx;   gy_f=gy;   gz_f=gz;
    lpf_init = true;
  } else {
    ax_f = lpf1(ax_f, ax_g, dt_s, ACC_LPF_HZ);
    ay_f = lpf1(ay_f, ay_g, dt_s, ACC_LPF_HZ);
    az_f = lpf1(az_f, az_g, dt_s, ACC_LPF_HZ);

    gx_f = lpf1(gx_f, gx, dt_s, GYRO_LPF_HZ);
    gy_f = lpf1(gy_f, gy, dt_s, GYRO_LPF_HZ);
    gz_f = lpf1(gz_f, gz, dt_s, GYRO_LPF_HZ);
  }

  // norms for logic
  float amag = sqrtf(ax_f*ax_f + ay_f*ay_f + az_f*az_f);
  float amag_err = fabsf(amag - 1.0f);
  float gyro_norm = sqrtf(gx_f*gx_f + gy_f*gy_f + gz_f*gz_f);

  // detect "still" (only for beta handling)
  static float still_timer = 0.0f;
  bool still_now = (amag_err < AMAG_STILL_ERR_G) && (gyro_norm < GYRO_STILL_RAD_S);
  if (still_now) still_timer += dt_s;
  else           still_timer  = 0.0f;

  bool still_flag = (still_timer >= STILL_HOLD_S);

  // detect rotation episode end -> settle boost
  static bool  rotating = false;
  static float settle_timer = 0.0f;

  if (gyro_norm > ROT_RATE_TRIG_RAD_S) {
    rotating = true;
    settle_timer = 0.0f;
  } else {
    if (rotating) {
      settle_timer += dt_s;
      if (settle_timer >= SETTLE_AFTER_ROT_S) {
        rotating = false;
        settle_timer = 0.0f;
      }
    }
  }

  bool settle_flag = rotating;

  // beta schedule
  float beta_used = BETA_DYNAMIC;

  if (amag_err > AMAG_THR_HARD) {
    beta_used = 0.0f; // gyro-only
  } else if (amag_err > AMAG_THR_SOFT) {
    beta_used = BETA_NOM * BETA_SOFT_SCALE;
  } else {
    if (still_flag || settle_flag) beta_used = BETA_SETTLE;
    else                           beta_used = BETA_DYNAMIC;
  }
  filter.beta = beta_used;

  // update Madgwick with filtered signals
  if (mx == 0.0f && my == 0.0f && mz == 0.0f) {
    filter.updateIMU(gx_f, gy_f, gz_f, ax_f, ay_f, az_f);
  } else {
    filter.update(gx_f, gy_f, gz_f, ax_f, ay_f, az_f, mx, my, mz);
  }

  // quaternion
  float qw = filter.q0, qx = filter.q1, qy = filter.q2, qz = filter.q3;

  // Euler
  float roll_rad, pitch_rad, yaw_rad;
  quatToEulerZYX(qw, qx, qy, qz, roll_rad, pitch_rad, yaw_rad);
  const float rad2deg = 180.0f / PI;
  float roll_deg  = roll_rad  * rad2deg;
  float pitch_deg = pitch_rad * rad2deg;
  float yaw_deg   = yaw_rad   * rad2deg;

  // gravity removal explicitly (BODY)
  float gxb, gyb, gzb;
  gravityBodyFromQuat(qw, qx, qy, qz, gxb, gyb, gzb);

  // linear accel in BODY (g)
  float ax_lin_b_g = ax_f - gxb;
  float ay_lin_b_g = ay_f - gyb;
  float az_lin_b_g = az_f - gzb;

  // to m/s^2
  const float g0 = 9.81f;
  float ax_lin_b = ax_lin_b_g * g0;
  float ay_lin_b = ay_lin_b_g * g0;
  float az_lin_b = az_lin_b_g * g0;

  // rotate to WORLD
  float ax_lin_w, ay_lin_w, az_lin_w;
  quatRotateBodyToWorld(qw, qx, qy, qz, ax_lin_b, ay_lin_b, az_lin_b, ax_lin_w, ay_lin_w, az_lin_w);

  uint32_t t_ms = millis();

  // ECO logic
  bool send_stream    = true;
  bool log_to_ram     = false;

  if (ECO_MODE == 1) {
    bool in_time_window = (t_ms >= START_LOG_MS) && (t_ms <= (START_LOG_MS + LOG_DURATION_MS));
    send_stream = in_time_window;
    log_to_ram  = in_time_window;
  } else if (ECO_MODE == 2) {
    bool in_time_window = (t_ms >= START_LOG_MS) && (t_ms <= (START_LOG_MS + LOG_DURATION_MS));
    send_stream = true;
    log_to_ram  = in_time_window;
  } else {
    send_stream = true;
    log_to_ram  = false;
  }

  // RAM log
  bool logging_now = log_to_ram && (g_logCount < LOG_MAX_SAMPLES);
  static bool was_logging = false;
  if (logging_now && !was_logging) Serial.println("=== START RAM LOG WINDOW ===");
  if (!logging_now && was_logging) Serial.println("=== END RAM LOG WINDOW ===");
  was_logging = logging_now;

  if (logging_now) {
    Sample &s = g_logBuf[g_logCount++];
    s.t_ms = t_ms;
    s.qw = qw; s.qx = qx; s.qy = qy; s.qz = qz;
    s.ax = ax_f; s.ay = ay_f; s.az = az_f;
    s.gx = gx_f; s.gy = gy_f; s.gz = gz_f;
    s.mx = mx; s.my = my; s.mz = mz;
  }

  if (!send_stream) return;

  // BLE payload
  float payload[11];
  payload[0]  = (float)t_ms;
  payload[1]  = qw;
  payload[2]  = qx;
  payload[3]  = qy;
  payload[4]  = qz;
  payload[5]  = ax_lin_w;
  payload[6]  = ay_lin_w;
  payload[7]  = az_lin_w;
  payload[8]  = gx_f;
  payload[9]  = gy_f;
  payload[10] = gz_f;
  telemetryUpdate(payload, 11);

  // --- update motion metric for NEXT cycle (drives dynamic cadence) ---
  float m_gyro = clamp01(gyro_norm / GYRO_NORM_FULL);
  float m_amg  = clamp01(amag_err  / AMAG_ERR_FULL);
  float motion_metric = (m_gyro > m_amg) ? m_gyro : m_amg;
  g_prevMotionMetric = motion_metric;

// ---------- Serial output (human-readable, organized) ----------
// Use this instead of the clean CSV block at the end of loop().

Serial.print("t_ms="); Serial.print(t_ms);

Serial.print(" | q=[");
Serial.print(qw, 4); Serial.print(", ");
Serial.print(qx, 4); Serial.print(", ");
Serial.print(qy, 4); Serial.print(", ");
Serial.print(qz, 4); Serial.print("]");

Serial.print(" | a_lin_w(m/s2)=[");
Serial.print(ax_lin_w, 2); Serial.print(", ");
Serial.print(ay_lin_w, 2); Serial.print(", ");
Serial.print(az_lin_w, 2); Serial.print("]");

Serial.print(" | w(rad/s)=[");
Serial.print(gx_f, 3); Serial.print(", ");
Serial.print(gy_f, 3); Serial.print(", ");
Serial.print(gz_f, 3); Serial.print("]");

Serial.print(" | euler(deg)=[");
Serial.print(roll_deg, 1); Serial.print(", ");
Serial.print(pitch_deg, 1); Serial.print(", ");
Serial.print(yaw_deg, 1); Serial.print("]");

Serial.print(" | amag=");     Serial.print(amag, 3);
Serial.print(" err=");        Serial.print(amag_err, 3);
Serial.print(" gyroN=");      Serial.print(gyro_norm, 3);

Serial.print(" | beta=");     Serial.print(beta_used, 3);
Serial.print(" still=");      Serial.print(still_flag ? 1 : 0);
Serial.print(" settle=");     Serial.print(settle_flag ? 1 : 0);

Serial.print(" | Hz tgt=");   Serial.print(g_targetFreqHz, 1);
Serial.print(" act=");        Serial.print(sampleFreq, 1);

Serial.println();

}
