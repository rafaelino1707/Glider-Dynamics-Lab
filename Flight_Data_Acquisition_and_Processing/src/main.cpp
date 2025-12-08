#include <Arduino.h>
#include <ArduinoBLE.h>
#include "imu.h"
#include "madgwick.h"
#include "orientation.h"
#include "telemetry.h"

// -----------------------------------------------------
// Modo de operação
// -----------------------------------------------------
// 1 = modo "eco": só envia entre 120 s e 150 s após o boot
// 0 = envia SEMPRE (sem janela de tempo)
#define ECO_MODE 1

// -----------------------------------------------------
// Configuração do Madgwick
// -----------------------------------------------------
MadgwickAHRS filter(100.0f, 0.1f); // 100 Hz, beta = 0.1

const float    SAMPLE_FREQ      = 100.0f;               // Hz
const uint32_t SAMPLE_PERIOD_US = 1000000UL / 100;      // 100 Hz

// -----------------------------------------------------
// Janela de logging (só usada se ECO_MODE == 1)
// -----------------------------------------------------
const uint32_t START_LOG_MS    = 120000UL;  // começa a logar aos 120 s
const uint32_t LOG_DURATION_MS = 30000UL;   // loga durante 30 s (até 150 s)

// -----------------------------------------------------
// Logging em RAM (backup se BLE falhar)
// Guarda: t_ms, q, accel bruta (g), gyro (rad/s), mag
// -----------------------------------------------------
struct Sample {
    uint32_t t_ms;
    float qw, qx, qy, qz;
    float ax, ay, az;      // em g (como vem da IMU)
    float gx, gy, gz;      // em rad/s
    float mx, my, mz;      // magnetómetro
};

// Ajusta se precisares (2000 ≈ 20 s a 100 Hz)
constexpr size_t LOG_MAX_SAMPLES = 2000;

Sample g_logBuf[LOG_MAX_SAMPLES];
size_t g_logCount = 0;

// Dump do log em RAM pela Serial (CSV simples)
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
    // --------- comando pela Serial: 'D' → dump RAM ---------
    if (Serial.available() > 0) {
        int c = Serial.read();
        if (c == 'D' || c == 'd') {
            Serial.println("=== RAM LOG DUMP BEGIN ===");
            dumpRamLogToSerial();
            Serial.println("=== RAM LOG DUMP END ===");
        }
    }

    // --------- manter BLE vivo SEMPRE ---------
    BLE.poll();

    // --------- controlo de frequência (100 Hz) ---------
    static uint32_t last = micros();
    uint32_t now = micros();
    uint32_t dt  = now - last;
    if (dt < SAMPLE_PERIOD_US) {
        // mantém o comportamento antigo
        delayMicroseconds(SAMPLE_PERIOD_US - dt);
        return;
    }
    last = micros();

    // --------- leitura IMU ---------
    float ax, ay, az;
    float gx, gy, gz;    // em rad/s (para Madgwick)
    float mx, my, mz;

    imuRead(ax, ay, az, gx, gy, gz, mx, my, mz);

    // --------- Atualizar Madgwick ---------
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

    // Gyro em deg/s (para debug / logging)
    float gx_deg = gx * rad2deg;
    float gy_deg = gy * rad2deg;
    float gz_deg = gz * rad2deg;

    // Aceleração linear (sem gravidade) em "g"
    float ax_lin, ay_lin, az_lin;
    removeGravityFromAccel(qw, qx, qy, qz, ax, ay, az, ax_lin, ay_lin, az_lin);

    // Conversão para m/s²
    const float g0 = 9.81f;

    float ax_mps2     = ax     * g0;
    float ay_mps2     = ay     * g0;
    float az_mps2     = az     * g0;
    float ax_lin_mps2 = ax_lin * g0;
    float ay_lin_mps2 = ay_lin * g0;
    float az_lin_mps2 = az_lin * g0;

    uint32_t t_ms = millis();

    // --------- Janela de logging / envio ---------
    bool in_logging_window = true;

#if ECO_MODE
    in_logging_window =
        (t_ms >= START_LOG_MS) &&
        (t_ms <= (START_LOG_MS + LOG_DURATION_MS));
#endif

    if (!in_logging_window) {
        // Fora da janela em modo ECO:
        // IMU + Madgwick correm e BLE.poll() é chamado,
        // mas não enviamos nada nem pela Serial nem por BLE.
        return;
    }

    // --------- logging em RAM (backup) ---------
    if (g_logCount < LOG_MAX_SAMPLES) {
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

    // --------- payload BLE completo (14 floats) ---------
    // [ t_ms, qw, qx, qy, qz, ax, ay, az, gx, gy, gz, mx, my, mz ]
    float payload[14];
    payload[0]  = (float)t_ms;  // timestamp em ms (float por simplicidade)
    payload[1]  = qw;
    payload[2]  = qx;
    payload[3]  = qy;
    payload[4]  = qz;
    payload[5]  = ax;           // accel em g (ou troca para ax_mps2 se preferires)
    payload[6]  = ay;
    payload[7]  = az;
    payload[8]  = gx;           // gyro em rad/s
    payload[9]  = gy;
    payload[10] = gz;
    payload[11] = mx;           // magnetómetro (µT, conforme devolve o imuRead)
    payload[12] = my;
    payload[13] = mz;

    telemetryUpdate(payload, 14);

    // --------- CSV na Serial ---------
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
