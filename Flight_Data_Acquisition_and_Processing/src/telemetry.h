#pragma once
#include <Arduino.h>
#include <ArduinoBLE.h>

bool telemetryBegin();
void telemetryUpdate(const float* data, size_t len); // len = número de floats
