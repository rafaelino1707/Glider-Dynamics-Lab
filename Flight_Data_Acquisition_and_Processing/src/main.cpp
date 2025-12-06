#include <Arduino.h>
#include <Arduino_BMI270_BMM150.h>  // IMU Rev2

void setup() {
  Serial.begin(115200);
  while (!Serial) {
    ; // esperar porta série
  }

  if (!IMU.begin()) {
    Serial.println("Failed to initialize IMU!");
    while (1);
  }

  Serial.print("Accel sample rate = ");
  Serial.print(IMU.accelerationSampleRate());
  Serial.println(" Hz");

  Serial.print("Gyro sample rate  = ");
  Serial.print(IMU.gyroscopeSampleRate());
  Serial.println(" Hz");
}

void loop() {
  float ax, ay, az;
  float gx, gy, gz;

  if (IMU.accelerationAvailable()) {
    IMU.readAcceleration(ax, ay, az);
  }

  if (IMU.gyroscopeAvailable()) {
    IMU.readGyroscope(gx, gy, gz);
  }

  Serial.print("a [m/s^2] = ");
  Serial.print(ax); Serial.print(", ");
  Serial.print(ay); Serial.print(", ");
  Serial.print(az);

  Serial.print("   g [deg/s] = ");
  Serial.print(gx); Serial.print(", ");
  Serial.print(gy); Serial.print(", ");
  Serial.println(gz);

  delay(10); // ~100 Hz
}
