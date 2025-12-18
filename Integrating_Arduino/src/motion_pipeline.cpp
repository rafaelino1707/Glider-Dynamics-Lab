#include "motion_pipeline.h"
#include "orientation.h"   // removeGravityFromAccel, rotateBodyToWorld, quatToEulerZYX
#include <math.h>

static inline float norm3(float x, float y, float z) {
  return sqrtf(x*x + y*y + z*z);
}

MotionPipeline::MotionPipeline() {
  reset();
}

void MotionPipeline::setConfig(const PipelineConfig& c) {
  cfg_ = c;

  acc_lpf_.setCutoff(cfg_.acc_lpf_hz);
  gyro_lpf_.setCutoff(cfg_.gyro_lpf_hz);

  acc_bias_w_.setTau(cfg_.acc_bias_tau_s);

  integ_.setConfig(cfg_.integ);
}

void MotionPipeline::reset() {
  acc_lpf_.reset({0,0,0});
  gyro_lpf_.reset({0,0,0});

  hx_.reset(0.0f);
  hy_.reset(0.0f);
  hz_.reset(0.0f);

  acc_bias_w_.reset({0,0,0});

  integ_.reset();

  still_timer_s_ = 0.0f;
  still_ = false;
}

bool MotionPipeline::computeStill(Vec3f acc_g, Vec3f gyro_rad, float dt_s) {
  const float wmag = norm3(gyro_rad.x, gyro_rad.y, gyro_rad.z);
  const float amag = norm3(acc_g.x, acc_g.y, acc_g.z);

  const bool cond =
      (wmag < cfg_.still.gyro_still_rad_s) &&
      (fabsf(amag - 1.0f) < cfg_.still.amag_still_err_g);

  if (cond) still_timer_s_ += dt_s;
  else      still_timer_s_  = 0.0f;

  return (still_timer_s_ >= cfg_.still.hold_s);
}

void MotionPipeline::update(float dt_s,
                            float qw, float qx, float qy, float qz,
                            Vec3f acc_g_body,
                            Vec3f gyro_rad_body,
                            PipelineOut& out)
{
  // dt protection
  if (!(dt_s > 0.0f) || !isfinite(dt_s)) dt_s = cfg_.integ.dt_min_s;
  dt_s = clampf(dt_s, cfg_.integ.dt_min_s, cfg_.integ.dt_max_s);

  // pre-LPF
  Vec3f acc_f = acc_lpf_.step(acc_g_body, dt_s);
  Vec3f gyr_f = gyro_lpf_.step(gyro_rad_body, dt_s);

  // Hampel on accel (spikes)
  acc_f.x = hx_.step(acc_f.x, cfg_.hampel_nsigma);
  acc_f.y = hy_.step(acc_f.y, cfg_.hampel_nsigma);
  acc_f.z = hz_.step(acc_f.z, cfg_.hampel_nsigma);

  // still detect (filtered)
  still_ = computeStill(acc_f, gyr_f, dt_s);

  // gravity removal -> a_lin em g (body)
  float ax_lin_g = 0.0f, ay_lin_g = 0.0f, az_lin_g = 0.0f;
  removeGravityFromAccel(qw, qx, qy, qz,
                         acc_f.x, acc_f.y, acc_f.z,
                         ax_lin_g, ay_lin_g, az_lin_g);

  // body->world (em g)
  float axw_g = 0.0f, ayw_g = 0.0f, azw_g = 0.0f;
  rotateBodyToWorld(qw, qx, qy, qz,
                    ax_lin_g, ay_lin_g, az_lin_g,
                    axw_g, ayw_g, azw_g);

  // converter para m/s^2
  Vec3f a_w = { axw_g * 9.80665f, ayw_g * 9.80665f, azw_g * 9.80665f };

  // bias estimator (apenas quando still) + remoção
  acc_bias_w_.update(a_w, dt_s, still_);
  Vec3f b = acc_bias_w_.bias();
  a_w.x -= b.x; a_w.y -= b.y; a_w.z -= b.z;

  // integração
  integ_.update(dt_s, a_w.x, a_w.y, a_w.z);

  // ZUPT hard (requer MotionIntegrator::setVelocity)
  bool zupt_applied = false;
  if (cfg_.still.enable_zupt && still_) {
    integ_.setVelocity(0.0f, 0.0f, 0.0f);
    zupt_applied = true;
  }

  // euler
  float roll = 0.0f, pitch = 0.0f, yaw = 0.0f;
  quatToEulerZYX(qw, qx, qy, qz, roll, pitch, yaw);

  // fill output
  out.qw = qw; out.qx = qx; out.qy = qy; out.qz = qz;
  out.roll_deg = roll; out.pitch_deg = pitch; out.yaw_deg = yaw;

  out.a_g_body = acc_f;
  out.w_rad_body = gyr_f;

  out.a_lin_g_body = { ax_lin_g, ay_lin_g, az_lin_g };
  out.a_lin_mps2_world = a_w;

  out.v_w = integ_.vel_w();
  out.p_w = integ_.pos_w();

  out.still = still_;
  out.zupt_applied = zupt_applied;
  out.rejected = integ_.lastRejected();
  out.jerk = integ_.lastJerk();
}
