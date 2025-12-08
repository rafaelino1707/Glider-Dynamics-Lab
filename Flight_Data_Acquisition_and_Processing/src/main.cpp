#include <Arduino.h>
#include <ArduinoBLE.h>
#include "imu.h"
#include "madgwick.h"
#include "orientation.h"
#include "telemetry.h"

// -----------------------------------------------------
// Operation Mode
// -----------------------------------------------------
// 1 = "eco" mode: only sends data between 120 s and 150 s since boot
// 0 = always sending data (no time window)
#define ECO_MODE 1

// -----------------------------------------------------
// Madgwick Configuration
// -----------------------------------------------------
MadgwickAHRS filter(100.0f, 0.1f); // 100 Hz, beta = 0.1

const float    SAMPLE_FREQ      = 100.0f;               // Hz
const uint32_t SAMPLE_PERIOD_US = 1000000UL / 100;      // 100 Hz

// -----------------------------------------------------
// Logging window (only used if ECO_MODE == 1)
// -----------------------------------------------------
const uint32_t START_LOG_MS    = 120000UL;  // Starts logging at 120 s
const uint32_t LOG_DURATION_MS = 30000UL;   // Logs for 30 s (until 150 s)

// -----------------------------------------------------
// RAM Logging (backup if BLE fails)
// Saves: t_ms, q, raw accel (g), gyro (rad/s), mag
// -----------------------------------------------------
struct Sample {
    uint32_t t_ms;
    float qw, qx, qy, qz;
    float ax, ay, az;      //  g 
    float gx, gy, gz;      //  rad/s
    float mx, my, mz;      // magnetometer
};

// RAM Log max samples 
constexpr size_t LOG_MAX_SAMPLES = 2000;

Sample g_logBuf[LOG_MAX_SAMPLES];
size_t g_logCount = 0;

// Dump RAM log to Serial Monitor (CSV)
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

void setup() {
    Serial.begin(115200);
    delay(2000);

    Serial.println("Glider IMU + Madgwick");

    if (!imuBegin()) {
        Serial.println("IMU init failed");
        while (1) {
            delay(1000);
        }
    }
    Serial.println("IMU OK");

    if (!telemetryBegin()) {
        Serial.println("BLE init failed");
        while (1) {
            delay(1000);
        }
    }
    Serial.println("BLE advertising");
}

void loop() {
    // --------- Serial Command: 'D' → dump RAM ---------
    if (Serial.available() > 0) {
        int c = Serial.read();
        if (c == 'D' || c == 'd') {
            Serial.println("=== RAM LOG DUMP BEGIN ===");
            dumpRamLogToSerial();
            Serial.println("=== RAM LOG DUMP END ===");
        }
    }

    // --------- Always keeping BLE alive ---------
    BLE.poll();

    // --------- Frequency control (100 Hz) ---------
    static uint32_t last = micros();
    uint32_t now = micros();
    uint32_t dt  = now - last;
    if (dt < SAMPLE_PERIOD_US) {
        delayMicroseconds(SAMPLE_PERIOD_US - dt);
        return;
    }
    last = micros();

    // --------- IMU Reading ---------
    float ax, ay, az;
    float gx, gy, gz;    // rad/s (for Madgwick)
    float mx, my, mz;

    imuRead(ax, ay, az, gx, gy, gz, mx, my, mz);

    // --------- Updating Madgwick ---------
    if (mx == 0.0f && my == 0.0f && mz == 0.0f) {
        filter.updateIMU(gx, gy, gz, ax, ay, az);
    } else {
        filter.update(gx, gy, gz, ax, ay, az, mx, my, mz);
    }

    // Actual Quaternion
    float qw = filter.q0;
    float qx = filter.q1;
    float qy = filter.q2;
    float qz = filter.q3;

    // Euler [rad]
    float roll_rad, pitch_rad, yaw_rad;
    quatToEulerZYX(qw, qx, qy, qz, roll_rad, pitch_rad, yaw_rad);

    // Euler [deg]
    const float rad2deg = 180.0f / PI;
    float roll_deg  = roll_rad  * rad2deg;
    float pitch_deg = pitch_rad * rad2deg;
    float yaw_deg   = yaw_rad   * rad2deg;

    // Gyro [deg/s] (debug / logging)
    float gx_deg = gx * rad2deg;
    float gy_deg = gy * rad2deg;
    float gz_deg = gz * rad2deg;

    // Linear Acceleration (without gravity)
    float ax_lin, ay_lin, az_lin;
    removeGravityFromAccel(qw, qx, qy, qz, ax, ay, az, ax_lin, ay_lin, az_lin);

    // Conversion to m/s²
    const float g0 = 9.81f;

    float ax_mps2     = ax     * g0;
    float ay_mps2     = ay     * g0;
    float az_mps2     = az     * g0;
    float ax_lin_mps2 = ax_lin * g0;
    float ay_lin_mps2 = ay_lin * g0;
    float az_lin_mps2 = az_lin * g0;

    uint32_t t_ms = millis();

    // --------- Logging / sending window ---------
    bool in_logging_window = true;

#if ECO_MODE
    // ECO_MODE = 1  → logs only between 120 s and 150 s
    in_logging_window =
        (t_ms >= START_LOG_MS) &&
        (t_ms <= (START_LOG_MS + LOG_DURATION_MS));
#endif

    // ----- RAM logger decision -----
    bool log_to_ram = false;

#if ECO_MODE
    // In flight (ECO_MODE=1) logs to RAM only inside the window
    log_to_ram = in_logging_window;
#else
    // In test mode (ECO_MODE=0) never uses RAM
    log_to_ram = false;
#endif

    // If outside logging window (in ECO_MODE=1),
    // do nothing (IMU + Madgwick still run).
    if (!in_logging_window) {
        return;
    }

    // --------- RAM logging (backup, only if log_to_ram=true) ---------
    if (log_to_ram && g_logCount < LOG_MAX_SAMPLES) {
        Sample &s = g_logBuf[g_logCount++];

        s.t_ms = t_ms;

        s.qw = qw;
        s.qx = qx;
        s.qy = qy;
        s.qz = qz;

        s.ax = ax;
        s.ay = ay;
        s.az = az;

        s.gx = gx;
        s.gy = gy;
        s.gz = gz;

        s.mx = mx;
        s.my = my;
        s.mz = mz;
    }

    // --------- Full BLE payload (14 floats) ---------
    // [ t_ms, qw, qx, qy, qz, ax, ay, az, gx, gy, gz, mx, my, mz ]
    float payload[14];
    payload[0]  = (float)t_ms;  // timestamp [ms]
    payload[1]  = qw;
    payload[2]  = qx;
    payload[3]  = qy;
    payload[4]  = qz;
    payload[5]  = ax;           // accel [g]
    payload[6]  = ay;
    payload[7]  = az;
    payload[8]  = gx;           // gyro [rad/s]
    payload[9]  = gy;
    payload[10] = gz;
    payload[11] = mx;           // magnetometer [µT]
    payload[12] = my;
    payload[13] = mz;

    telemetryUpdate(payload, 14);

    // --------- CSV over Serial ---------
    Serial.print(t_ms);        Serial.print(',');
    Serial.print(qw, 6);       Serial.print(',');
    Serial.print(qx, 6);       Serial.print(',');
    Serial.print(qy, 6);       Serial.print(',');
    Serial.print(qz, 6);       Serial.print(',');

    Serial.print(ax_mps2, 6);  Serial.print(',');
    Serial.print(ay_mps2, 6);  Serial.print(',');
    Serial.print(az_mps2, 6);  Serial.print(',');

    Serial.print(gx_deg, 6);   Serial.print(',');
    Serial.print(gy_deg, 6);   Serial.print(',');
    Serial.print(gz_deg, 6);   Serial.print(',');

    Serial.print(roll_deg, 6);   Serial.print(',');
    Serial.print(pitch_deg, 6);  Serial.print(',');
    Serial.print(yaw_deg, 6);    Serial.print(',');

    Serial.print(ax_lin_mps2, 6); Serial.print(',');
    Serial.print(ay_lin_mps2, 6); Serial.print(',');
    Serial.print(az_lin_mps2, 6);

    Serial.println();
}
