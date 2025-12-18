#include "telemetry.h"
#include <ArduinoBLE.h>
#include <string.h>

// -----------------------------------------------------
// Connection State BLE
// -----------------------------------------------------
static bool g_bleConnected = false;

// UUIDs custom
static BLEService imuService("12345678-1234-1234-1234-1234567890AB");

// Char 1: IMU legacy (14 floats = 56 bytes)
static BLECharacteristic imuChar(
    "12345678-1234-1234-1234-1234567890AC",
    BLERead | BLENotify,
    56
);

// Char 2: NAV (10 floats = 40 bytes)
static BLECharacteristic navChar(
    "12345678-1234-1234-1234-1234567890AD",
    BLERead | BLENotify,
    40
);

// callbacks prototypes
static void onBleConnect(BLEDevice central);
static void onBleDisconnect(BLEDevice central);

// -----------------------------------------------------
// Telemetry/BLE initialization
// -----------------------------------------------------
bool telemetryBegin() {
    if (!BLE.begin()) {
        return false;
    }

    BLE.setLocalName("GliderIMU");
    BLE.setDeviceName("GliderIMU");
    BLE.setAdvertisedService(imuService);

    imuService.addCharacteristic(imuChar);
    imuService.addCharacteristic(navChar);
    BLE.addService(imuService);

    // Initial values
    uint8_t zero56[56] = {0};
    uint8_t zero40[40] = {0};
    imuChar.writeValue(zero56, sizeof(zero56));
    navChar.writeValue(zero40, sizeof(zero40));

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
void telemetryUpdateIMU14(const float* data14) {
    if (!data14) return;
    uint8_t buf[56];
    memcpy(buf, data14, 56);
    imuChar.writeValue(buf, 56);
    BLE.poll();
}

void telemetryUpdateNAV10(const float* data10) {
    if (!data10) return;
    uint8_t buf[40];
    memcpy(buf, data10, 40);
    navChar.writeValue(buf, 40);
    BLE.poll();
}

// -----------------------------------------------------
// Connection State BLE
// -----------------------------------------------------
bool telemetryIsConnected() {
    return g_bleConnected;
}

// -----------------------------------------------------
// Callbacks
// -----------------------------------------------------
static void onBleConnect(BLEDevice central) {
    (void)central;
    g_bleConnected = true;
}

static void onBleDisconnect(BLEDevice central) {
    (void)central;
    g_bleConnected = false;
}

void telemetryPoll() {
    BLE.poll();
}
