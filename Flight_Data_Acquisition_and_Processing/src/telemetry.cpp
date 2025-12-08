#include "telemetry.h"
#include <ArduinoBLE.h>

// -----------------------------------------------------
// Estado da ligação BLE
// -----------------------------------------------------
static bool g_bleConnected = false;

// UUIDs custom
static BLEService imuService("12345678-1234-1234-1234-1234567890AB");

static BLECharacteristic imuChar(
    "12345678-1234-1234-1234-1234567890AC",
    BLERead | BLENotify,
    56    // 14 floats * 4 bytes
);

// Protótipos dos callbacks
static void onBleConnect(BLEDevice central);
static void onBleDisconnect(BLEDevice central);

// -----------------------------------------------------
// Inicialização do BLE / Telemetria
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

    // valor inicial (zeros)
    uint8_t zero[56] = {0};
    imuChar.writeValue(zero, sizeof(zero));

    // callbacks de ligação
    BLE.setEventHandler(BLEConnected,    onBleConnect);
    BLE.setEventHandler(BLEDisconnected, onBleDisconnect);

    g_bleConnected = false;
    BLE.advertise();
    return true;
}

// -----------------------------------------------------
// Enviar amostra por BLE
// -----------------------------------------------------
void telemetryUpdate(const float* data, size_t len) {
    if (len == 0) return;

    const size_t bytes = len * sizeof(float);
    if (bytes > 56) return; // 14 floats no máximo

    uint8_t buf[56];
    memcpy(buf, data, bytes);

    imuChar.writeValue(buf, bytes);
    BLE.poll();
}

// -----------------------------------------------------
// Estado da ligação BLE (usado no loop principal)
// -----------------------------------------------------
bool telemetryIsConnected() {
    return g_bleConnected;
}

// -----------------------------------------------------
// Callbacks de ligação/desligação
// -----------------------------------------------------
static void onBleConnect(BLEDevice central) {
    g_bleConnected = true;
}

static void onBleDisconnect(BLEDevice central) {
    g_bleConnected = false;
}
