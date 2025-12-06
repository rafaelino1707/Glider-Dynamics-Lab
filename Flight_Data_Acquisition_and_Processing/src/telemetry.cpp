#include "telemetry.h"

static BLEService imuService("180D"); // UUID qualquer
// 7 floats = 28 bytes: qw,qx,qy,qz,ax,ay,az
static BLECharacteristic imuChar("2A37",
                                 BLERead | BLENotify,
                                 28);

bool telemetryBegin() {
    if (!BLE.begin()) {
        return false;
    }

    BLE.setLocalName("GliderIMU");
    BLE.setAdvertisedService(imuService);

    imuService.addCharacteristic(imuChar);
    BLE.addService(imuService);

    BLE.advertise();
    return true;
}

void telemetryUpdate(const float* data, size_t len) {
    if (!BLE.connected()) return;
    if (len == 0) return;

    const size_t bytes = len * sizeof(float);
    if (bytes > 28) return; // proteger

    uint8_t buf[28];
    memcpy(buf, data, bytes);
    imuChar.writeValue(buf, bytes);
}
