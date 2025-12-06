#include <Arduino.h>
#include "imu.h"
#include "madgwick.h"
#include "orientation.h"
#include "telemetry.h"   // se já o tens; se não, remove esta linha e as partes de BLE

MadgwickAHRS filter(100.0f, 0.1f); // 100 Hz, beta = 0.1

const float SAMPLE_FREQ = 100.0f;                     // Hz
const uint32_t SAMPLE_PERIOD_US = 1000000UL / 100;   // 100 Hz

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

    // Se estiveres a usar BLE:
    if (!telemetryBegin()) {
        Serial.println("BLE init failed");
        while (1) {
            delay(1000);
        }
    }
    Serial.println("BLE advertising");
}

void loop() {
    static uint32_t last = micros();
    uint32_t now = micros();
    uint32_t dt = now - last;
    if (dt < SAMPLE_PERIOD_US) {
        delayMicroseconds(SAMPLE_PERIOD_US - dt);
        return;
    }
    last = micros();

    float ax, ay, az;
    float gx, gy, gz;    // em rad/s (para Madgwick)
    float mx, my, mz;

    imuRead(ax, ay, az, gx, gy, gz, mx, my, mz);

    // Atualizar Madgwick
    if (mx == 0.0f && my == 0.0f && mz == 0.0f) {
        filter.updateIMU(gx, gy, gz, ax, ay, az);
    } else {
        filter.update(gx, gy, gz, ax, ay, az, mx, my, mz);
    }

    // Quaternion atual
    float qw = filter.q0;
    float qx = filter.q1;
    float qy = filter.q2;
    float qz = filter.q3;

    // Euler (rad)
    float roll_rad, pitch_rad, yaw_rad;
    quatToEulerZYX(qw, qx, qy, qz, roll_rad, pitch_rad, yaw_rad);

    // Euler em graus
    const float rad2deg = 180.0f / PI;
    float roll_deg  = roll_rad  * rad2deg;
    float pitch_deg = pitch_rad * rad2deg;
    float yaw_deg   = yaw_rad   * rad2deg;

    // Gyro em deg/s para impressão (internamente está em rad/s)
    float gx_deg = gx * rad2deg;
    float gy_deg = gy * rad2deg;
    float gz_deg = gz * rad2deg;

    // Aceleração linear (sem gravidade) em "g"
    float ax_lin, ay_lin, az_lin;
    removeGravityFromAccel(qw, qx, qy, qz, ax, ay, az, ax_lin, ay_lin, az_lin);

    // Pacote simples para BLE (mantendo 7 floats: q + accel bruta)
    float payload[7];
    payload[0] = qw;
    payload[1] = qx;
    payload[2] = qy;
    payload[3] = qz;
    payload[4] = ax;
    payload[5] = ay;
    payload[6] = az;
    telemetryUpdate(payload, 7);

    // Logging pela porta série
    Serial.print("q = ");
    Serial.print(qw, 4); Serial.print(", ");
    Serial.print(qx, 4); Serial.print(", ");
    Serial.print(qy, 4); Serial.print(", ");
    Serial.print(qz, 4);

    Serial.print("   a = ");
    Serial.print(ax, 3); Serial.print(", ");
    Serial.print(ay, 3); Serial.print(", ");
    Serial.print(az, 3);

    Serial.print("   g[deg/s] = ");
    Serial.print(gx_deg, 2); Serial.print(", ");
    Serial.print(gy_deg, 2); Serial.print(", ");
    Serial.print(gz_deg, 2);

    Serial.print("   rpy[deg] = ");
    Serial.print(roll_deg, 2); Serial.print(", ");
    Serial.print(pitch_deg, 2); Serial.print(", ");
    Serial.print(yaw_deg, 2);

    Serial.print("   a_lin = ");
    Serial.print(ax_lin, 3); Serial.print(", ");
    Serial.print(ay_lin, 3); Serial.print(", ");
    Serial.print(az_lin, 3);

    Serial.println();
}
