#pragma once
#include <Arduino.h>
#include <ArduinoBLE.h>

bool telemetryBegin();
bool telemetryIsConnected();

// Mantém o pacote antigo (14 floats = 56 bytes)
void telemetryUpdateIMU14(const float* data14); // len fixo 14

// Novo pacote NAV (10 floats = 40 bytes): v/p + flags + jerk
void telemetryUpdateNAV10(const float* data10); // len fixo 10

void telemetryPoll();
