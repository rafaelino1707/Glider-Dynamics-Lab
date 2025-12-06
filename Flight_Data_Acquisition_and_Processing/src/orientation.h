#pragma once
#include <Arduino.h>

// Converte quaternion para ângulos de Euler (convenção ZYX: yaw, pitch, roll)
void quatToEulerZYX(float q0, float q1, float q2, float q3,
                    float& roll, float& pitch, float& yaw);

// Remove componente da gravidade da aceleração medida
// Assume que quando parado o módulo de 'a' é ≈ 1 (unidades "g" como estás a ver).
void removeGravityFromAccel(float q0, float q1, float q2, float q3,
                            float ax, float ay, float az,
                            float& ax_lin, float& ay_lin, float& az_lin);
