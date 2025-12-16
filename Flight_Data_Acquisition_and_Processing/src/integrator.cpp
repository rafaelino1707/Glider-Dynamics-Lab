#include "integrator.h"

MotionIntegrator::MotionIntegrator() {}

void MotionIntegrator::setConfig(const IntegratorConfig& cfg) {
  cfg_ = cfg;
}

void MotionIntegrator::reset() {
  x_prev_ = {0,0,0};
  y_prev_ = {0,0,0};
  a_out_  = {0,0,0};
  a_prev_ = {0,0,0};
  v_prev_ = {0,0,0};
  v_      = {0,0,0};
  p_      = {0,0,0};
  inited_ = false;

  last_rejected_ = false;
  last_jerk_ = 0.0f;
}

void MotionIntegrator::update(float dt_s, float ax_w, float ay_w, float az_w) {
  last_rejected_ = false;

  // dt clamp
  if (dt_s < cfg_.dt_min_s) dt_s = cfg_.dt_min_s;
  if (dt_s > cfg_.dt_max_s) dt_s = cfg_.dt_max_s;

  // first call init: set prevs and do nothing aggressive
  if (!inited_) {
    x_prev_ = {ax_w, ay_w, az_w};
    y_prev_ = {0,0,0}; // HPF starts at 0
    a_out_  = {0,0,0};
    a_prev_ = a_out_;
    v_prev_ = v_;
    inited_ = true;
    last_jerk_ = 0.0f;
    return;
  }

  // --- jerk gate (computed on RAW accel in world) ---
  float dax = (ax_w - x_prev_.x) / dt_s;
  float day = (ay_w - x_prev_.y) / dt_s;
  float daz = (az_w - x_prev_.z) / dt_s;
  float jerk = sqrtf(dax*dax + day*day + daz*daz);
  last_jerk_ = jerk;

  if (cfg_.jerk_max_mps3 > 0.0f && jerk > cfg_.jerk_max_mps3) {
    // Reject this sample: update prev input so we don't keep “jerk” forever,
    // but DO NOT integrate.
    x_prev_ = {ax_w, ay_w, az_w};
    last_rejected_ = true;
    return;
  }

  // --- HPF on accel (offset/drift killer) ---
  float ax_h = hpf1_step(y_prev_.x, ax_w, x_prev_.x, dt_s, cfg_.acc_hpf_hz);
  float ay_h = hpf1_step(y_prev_.y, ay_w, x_prev_.y, dt_s, cfg_.acc_hpf_hz);
  float az_h = hpf1_step(y_prev_.z, az_w, x_prev_.z, dt_s, cfg_.acc_hpf_hz);

  y_prev_ = {ax_h, ay_h, az_h};
  x_prev_ = {ax_w, ay_w, az_w};

  a_out_ = {ax_h, ay_h, az_h};

  // --- integrate accel -> vel (trapezoidal) ---
  v_prev_ = v_;
  v_.x += 0.5f * (a_prev_.x + a_out_.x) * dt_s;
  v_.y += 0.5f * (a_prev_.y + a_out_.y) * dt_s;
  v_.z += 0.5f * (a_prev_.z + a_out_.z) * dt_s;

  // --- velocity leak (drift killer) ---
  float leak = expLeak(dt_s, cfg_.vel_leak_tau_s);
  v_.x *= leak;
  v_.y *= leak;
  v_.z *= leak;

  // --- integrate vel -> pos (trapezoidal) ---
  p_.x += 0.5f * (v_prev_.x + v_.x) * dt_s;
  p_.y += 0.5f * (v_prev_.y + v_.y) * dt_s;
  p_.z += 0.5f * (v_prev_.z + v_.z) * dt_s;

  a_prev_ = a_out_;

  clampStates();
}

void MotionIntegrator::clampStates() {
  // optional safety clamps
  float vmax = cfg_.v_abs_max_mps;
  if (vmax > 0.0f) {
    v_.x = clampf(v_.x, -vmax, vmax);
    v_.y = clampf(v_.y, -vmax, vmax);
    v_.z = clampf(v_.z, -vmax, vmax);
  }
  float pmax = cfg_.p_abs_max_m;
  if (pmax > 0.0f) {
    p_.x = clampf(p_.x, -pmax, pmax);
    p_.y = clampf(p_.y, -pmax, pmax);
    p_.z = clampf(p_.z, -pmax, pmax);
  }
}
