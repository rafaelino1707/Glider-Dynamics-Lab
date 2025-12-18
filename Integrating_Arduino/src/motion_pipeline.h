#pragma once
#include <Arduino.h>
#include "integrator.h"
#include "robust_filters.h"

struct StillConfig {
  // still detect (raw)
  float gyro_still_rad_s = 0.08f;     // |w| < ...
  float amag_still_err_g = 0.05f;     // | |a|-1 | < ...
  float hold_s           = 0.35f;     // tem de aguentar isto

  // ZUPT
  bool  enable_zupt      = true;
};

struct PipelineConfig {
  // prefilters
  float acc_lpf_hz  = 12.0f;
  float gyro_lpf_hz = 12.0f;

  // hampel
  float hampel_nsigma = 4.0f;

  // bias estimator (no accel linear em world) quando still
  float acc_bias_tau_s = 20.0f;

  StillConfig still;
  IntegratorConfig integ;
};

struct PipelineOut {
  // orientation
  float qw,qx,qy,qz;
  float roll_deg,pitch_deg,yaw_deg;

  // sensors
  Vec3f a_g_body;
  Vec3f w_rad_body;

  // linear accel
  Vec3f a_lin_g_body;      // g
  Vec3f a_lin_mps2_world;  // m/s^2 (já com bias removida + filtros)

  // states
  Vec3  v_w;
  Vec3  p_w;

  // debug
  bool  still;
  bool  zupt_applied;
  bool  rejected;
  float jerk;
};

class MotionPipeline {
public:
  MotionPipeline();

  void setConfig(const PipelineConfig& c);
  void reset();

  // inputs:
  //  - quaternion (Madgwick)
  //  - accel em g (body)
  //  - gyro em rad/s (body)
  //  - dt em s
  void update(float dt_s,
              float qw,float qx,float qy,float qz,
              Vec3f acc_g_body,
              Vec3f gyro_rad_body,
              PipelineOut& out);

private:
  PipelineConfig cfg_;

  Vec3LPF acc_lpf_, gyro_lpf_;

  Hampel<9> hx_, hy_, hz_; // accel spikes por eixo

  BiasEstimator3 acc_bias_w_;

  MotionIntegrator integ_;

  // still hold
  float still_timer_s_ = 0.0f;
  bool  still_ = false;

  bool computeStill(Vec3f acc_g, Vec3f gyro_rad, float dt_s);
};
