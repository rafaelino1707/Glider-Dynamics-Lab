#include "imu.h"
#include <Arduino_BMI270_BMM150.h>

// Wrapper simples

bool imuBegin() {
    if (!IMU.begin()) {
        return false;
    }
    return true;
}

void imuRead(float& ax, float& ay, float& az,
             float& gx, float& gy, float& gz,
             float& mx, float& my, float& mz) {

    // initialize with zeros
    ax = ay = az = 0.0f;
    gx = gy = gz = 0.0f;
    mx = my = mz = 0.0f;

    if (IMU.accelerationAvailable()) {
        IMU.readAcceleration(ax, ay, az);
        // BMI270 library already in m/s^2
    }

    if (IMU.gyroscopeAvailable()) {
        IMU.readGyroscope(gx, gy, gz);
        // Library gives deg/s -> convert to rad/s
        const float deg2rad = PI / 180.0f;
        gx *= deg2rad;
        gy *= deg2rad;
        gz *= deg2rad;
    }

    if (IMU.magneticFieldAvailable()) {
        IMU.readMagneticField(mx, my, mz);
        // Em microTesla, ok para Madgwick
    }
}
