#pragma once
#include <Arduino.h>

bool imuBegin();
void imuRead(float& ax, float& ay, float& az,
             float& gx, float& gy, float& gz,
             float& mx, float& my, float& mz);
