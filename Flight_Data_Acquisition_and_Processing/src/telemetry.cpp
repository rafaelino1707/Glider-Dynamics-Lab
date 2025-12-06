#include "telemetry.h"
#include <ArduinoBLE.h>

// UUIDs custom
static BLEService imuService("12345678-1234-1234-1234-1234567890AB");

static BLECharacteristic imuChar(
    "12345678-1234-1234-1234-1234567890AC",
    BLERead | BLENotify,
    28    // 7 floats * 4 bytes
);

bool telemetryBegin() {
    if (!BLE.begin()) {
        return false;
    }

    BLE.setLocalName("GliderIMU");
    BLE.setDeviceName("GliderIMU");
    BLE.setAdvertisedService(imuService);

    imuService.addCharacteristic(imuChar);
    BLE.addService(imuService);

    // valor inicial
    uint8_t zero[4] = {0};
    imuChar.writeValue(zero, 4);

    BLE.advertise();
    return true;
}

void telemetryUpdate(const float* data, size_t len) {
    if (len == 0) return;

    const size_t bytes = len * sizeof(float);
    if (bytes > 28) return; // 7 floats no máximo

    uint8_t buf[28];
    memcpy(buf, data, bytes);

    imuChar.writeValue(buf, bytes);
    BLE.poll();  // ajuda a manter a stack BLE viva
}
