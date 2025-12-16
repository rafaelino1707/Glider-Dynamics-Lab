#pragma once
#include <Arduino.h>

struct Vec3 {
  float x, y, z;
};

struct IntegratorConfig {
  // --- dt safety ---
  float dt_min_s = 1e-4f;
  float dt_max_s = 0.20f;

  // --- accel HPF (offset/drift killer) ---
  // 0.05..0.20 Hz típico (quanto maior, mais “mata” drift mas também mata low-speed motion)
  float acc_hpf_hz = 0.10f;

  // --- velocity leak (drift killer sem ZUPT) ---
  // v = v * exp(-dt/tau). tau pequeno -> mata drift mais agressivo
  float vel_leak_tau_s = 6.0f;

  // --- jerk gate ---
  // se |da/dt| > jerk_max -> ignora update (ou aplica fallback)
  float jerk_max_mps3 = 80.0f;

  // --- optional: hard clamp de velocidade/pos (proteção) ---
  float v_abs_max_mps = 150.0f;   // grande, só para safety
  float p_abs_max_m   = 1e6f;     // safety
};

class MotionIntegrator {
public:
  MotionIntegrator();

  void setConfig(const IntegratorConfig& cfg);
  void reset();

  // Update with linear accel in WORLD frame (m/s^2) and dt (s)
  void update(float dt_s, float ax_w, float ay_w, float az_w);

  // Accessors
  Vec3 acc_w() const { return a_out_; }   // pós-HPF (o que foi integrado)
  Vec3 vel_w() const { return v_; }
  Vec3 pos_w() const { return p_; }

  // Debug
  bool lastRejected() const { return last_rejected_; }
  float lastJerk() const { return last_jerk_; }

private:
  IntegratorConfig cfg_;

  // HPF state per axis: y = HPF(x)
  Vec3 x_prev_{0,0,0};
  Vec3 y_prev_{0,0,0};

  // output accel after HPF
  Vec3 a_out_{0,0,0};

  // trapezoid integration memory
  Vec3 a_prev_{0,0,0};  // last a_out
  Vec3 v_prev_{0,0,0};  // last v

  // integrated states
  Vec3 v_{0,0,0};
  Vec3 p_{0,0,0};

  bool inited_ = false;

  // debug
  bool  last_rejected_ = false;
  float last_jerk_ = 0.0f;

  // helpers
  static inline float clampf(float x, float lo, float hi) {
    return (x < lo) ? lo : (x > hi ? hi : x);
  }

  // 1st-order HPF (discrete) with cutoff fc:
  // y[n] = a*(y[n-1] + x[n] - x[n-1]), with a = RC/(RC+dt), RC = 1/(2*pi*fc)
  static inline float hpf1_step(float y_prev, float x, float x_prev, float dt, float fc_hz) {
    if (fc_hz <= 0.0f || dt <= 0.0f) return x; // “no HPF”
    float RC = 1.0f / (2.0f * PI * fc_hz);
    float a  = RC / (RC + dt);
    return a * (y_prev + x - x_prev);
  }

  static inline float expLeak(float dt, float tau) {
    if (tau <= 0.0f) return 1.0f;
    return expf(-dt / tau);
  }

  void clampStates();
};
