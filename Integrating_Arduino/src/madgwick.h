#pragma once
#include <Arduino.h>

class MadgwickAHRS {
public:
    float beta;        // algorithm gain
    float sampleFreq;  // sample frequency in Hz

    float q0, q1, q2, q3; // quaternion

    MadgwickAHRS(float sampleFrequency = 100.0f, float betaDef = 0.1f);

    // Full 9-DOF update
    void update(float gx, float gy, float gz,
                float ax, float ay, float az,
                float mx, float my, float mz);

    // IMU-only update (sem magnetómetro)
    void updateIMU(float gx, float gy, float gz,
                   float ax, float ay, float az);
};
