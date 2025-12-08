#include "telemetry.h"
#include <ArduinoBLE.h>

// -----------------------------------------------------
// Connection State BLE
// -----------------------------------------------------
static bool g_bleConnected = false;

// UUIDs custom
static BLEService imuService("12345678-1234-1234-1234-1234567890AB");

static BLECharacteristic imuChar(
    "12345678-1234-1234-1234-1234567890AC",
    BLERead | BLENotify,
    56    // 14 floats * 4 bytes
);

// callbacks prototypes
static void onBleConnect(BLEDevice central);
static void onBleDisconnect(BLEDevice central);

// -----------------------------------------------------
// Telemetry/BLE inicialization
// -----------------------------------------------------
bool telemetryBegin() {
    if (!BLE.begin()) {
        return false;
    }

    BLE.setLocalName("GliderIMU");
    BLE.setDeviceName("GliderIMU");
    BLE.setAdvertisedService(imuService);

    imuService.addCharacteristic(imuChar);
    BLE.addService(imuService);

    // Initial values
    uint8_t zero[56] = {0};
    imuChar.writeValue(zero, sizeof(zero));

    // Connection callbacks
    BLE.setEventHandler(BLEConnected,    onBleConnect);
    BLE.setEventHandler(BLEDisconnected, onBleDisconnect);

    g_bleConnected = false;
    BLE.advertise();
    return true;
}

// -----------------------------------------------------
// Sending data through BLE
// -----------------------------------------------------
void telemetryUpdate(const float* data, size_t len) {
    if (len == 0) return;

    const size_t bytes = len * sizeof(float);
    if (bytes > 56) return; // Maximum of 14 floats

    uint8_t buf[56];
    memcpy(buf, data, bytes);

    imuChar.writeValue(buf, bytes);
    BLE.poll();
}

// -----------------------------------------------------
// Connectin State BLE 
// -----------------------------------------------------
bool telemetryIsConnected() {
    return g_bleConnected;
}

// -----------------------------------------------------
// Callbacks of Connection
// -----------------------------------------------------
static void onBleConnect(BLEDevice central) {
    g_bleConnected = true;
}

static void onBleDisconnect(BLEDevice central) {
    g_bleConnected = false;
}
