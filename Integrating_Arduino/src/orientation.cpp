#include "orientation.h"

void quatToEulerZYX(float q0, float q1, float q2, float q3,
                    float& roll, float& pitch, float& yaw)
{
    // Convenção ZYX (yaw-pitch-roll)
    // roll (x-axis rotation)
    float sinr_cosp = 2.0f * (q0 * q1 + q2 * q3);
    float cosr_cosp = 1.0f - 2.0f * (q1 * q1 + q2 * q2);
    roll = atan2f(sinr_cosp, cosr_cosp);

    // pitch (y-axis rotation)
    float sinp = 2.0f * (q0 * q2 - q3 * q1);
    if (fabsf(sinp) >= 1.0f) {
        // saturar a ±90 deg
        pitch = copysignf(PI / 2.0f, sinp);
    } else {
        pitch = asinf(sinp);
    }

    // yaw (z-axis rotation)
    float siny_cosp = 2.0f * (q0 * q3 + q1 * q2);
    float cosy_cosp = 1.0f - 2.0f * (q2 * q2 + q3 * q3);
    yaw = atan2f(siny_cosp, cosy_cosp);
}

// Gravidade prevista no referencial do sensor, a partir do quaternion
// (estas são as mesmas expressões que o Madgwick usa internamente)
static void gravityBodyFromQuat(float q0, float q1, float q2, float q3,
                                float& gx, float& gy, float& gz)
{
    gx = 2.0f * (q1 * q3 - q0 * q2);
    gy = 2.0f * (q0 * q1 + q2 * q3);
    gz = q0*q0 - q1*q1 - q2*q2 + q3*q3;
}

void removeGravityFromAccel(float q0, float q1, float q2, float q3,
                            float ax, float ay, float az,
                            float& ax_lin, float& ay_lin, float& az_lin)
{
    // Direção da gravidade no referencial do sensor (módulo ≈ 1)
    float gx, gy, gz;
    gravityBodyFromQuat(q0, q1, q2, q3, gx, gy, gz);

    // Assume-se que as unidades do acelerómetro são "g"
    // (como tens agora: módulo ≈ 1 parado).
    // Logo, aceleração linear ≈ medido - gravidade
    ax_lin = ax - gx;
    ay_lin = ay - gy;
    az_lin = az - gz;

}

void rotateBodyToWorld(float qw, float qx, float qy, float qz,
                              float bx, float by, float bz,
                              float &wx, float &wy, float &wz)
{
    // v_world = q ⊗ v_body ⊗ q_conj
    const float tx = 2.0f * (qy*bz - qz*by);
    const float ty = 2.0f * (qz*bx - qx*bz);
    const float tz = 2.0f * (qx*by - qy*bx);

    wx = bx + qw*tx + (qy*tz - qz*ty);
    wy = by + qw*ty + (qz*tx - qx*tz);
    wz = bz + qw*tz + (qx*ty - qy*tx);
}
