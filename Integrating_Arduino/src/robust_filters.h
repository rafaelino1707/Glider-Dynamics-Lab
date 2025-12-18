#pragma once
#include <Arduino.h>

struct Vec3f { float x, y, z; };

static inline float clampf(float v, float lo, float hi) { return (v < lo) ? lo : (v > hi ? hi : v); }

class OnePoleLPF {
public:
  OnePoleLPF() = default;
  void setCutoff(float fc_hz) { fc_ = fc_hz; }
  void reset(float x0=0.0f) { y_ = x0; inited_ = false; }
  float step(float x, float dt_s) {
    if (fc_ <= 0.0f || dt_s <= 0.0f) { y_ = x; inited_ = true; return y_; }
    if (!inited_) { y_ = x; inited_ = true; return y_; }
    const float w = 2.0f * PI * fc_;
    const float a = (w * dt_s) / (1.0f + w * dt_s);
    y_ += a * (x - y_);
    return y_;
  }
private:
  float fc_ = 0.0f;
  float y_ = 0.0f;
  bool  inited_ = false;
};

class Vec3LPF {
public:
  void setCutoff(float fc_hz){ fx_.setCutoff(fc_hz); fy_.setCutoff(fc_hz); fz_.setCutoff(fc_hz); }
  void reset(Vec3f v={0,0,0}){ fx_.reset(v.x); fy_.reset(v.y); fz_.reset(v.z); }
  Vec3f step(Vec3f v, float dt_s){
    return { fx_.step(v.x, dt_s), fy_.step(v.y, dt_s), fz_.step(v.z, dt_s) };
  }
private:
  OnePoleLPF fx_, fy_, fz_;
};

// Hampel filter (median + MAD) com janela pequena (N ímpar).
template<int N>
class Hampel {
  static_assert((N % 2) == 1, "N must be odd");
public:
  void reset(float x0=0.0f){
    for (int i=0;i<N;i++) buf_[i]=x0;
    idx_=0; filled_=true;
  }
  // nsigma típico 3..5
  float step(float x, float nsigma=4.0f){
    buf_[idx_] = x;
    idx_ = (idx_ + 1) % N;
    if (!filled_ && idx_==0) filled_=true;

    // copiar e ordenar (insertion sort) -> N pequeno
    float tmp[N];
    for (int i=0;i<N;i++) tmp[i]=buf_[i];
    for (int i=1;i<N;i++){
      float key=tmp[i];
      int j=i-1;
      while (j>=0 && tmp[j]>key){ tmp[j+1]=tmp[j]; j--; }
      tmp[j+1]=key;
    }
    float med = tmp[N/2];

    // MAD
    float dev[N];
    for (int i=0;i<N;i++) dev[i]=fabsf(buf_[i]-med);
    for (int i=1;i<N;i++){
      float key=dev[i];
      int j=i-1;
      while (j>=0 && dev[j]>key){ dev[j+1]=dev[j]; j--; }
      dev[j+1]=key;
    }
    float mad = dev[N/2];
    float sigma = 1.4826f * mad + 1e-6f;

    float z = fabsf(x - med) / sigma;
    return (z > nsigma) ? med : x;
  }
private:
  float buf_[N];
  int idx_ = 0;
  bool filled_ = false;
};

// Bias estimator (só atualiza quando still==true)
class BiasEstimator3 {
public:
  void setTau(float tau_s){ tau_ = tau_s; }
  void reset(Vec3f b={0,0,0}){ b_=b; inited_=false; }
  Vec3f bias() const { return b_; }

  // quando still: b <- LPF(meas)  (tau grande = lento)
  void update(Vec3f meas, float dt_s, bool still){
    if (!still || dt_s <= 0.0f) return;
    if (!inited_) { b_ = meas; inited_ = true; return; }
    float a = 1.0f - expf(-dt_s / max(1e-3f, tau_));
    b_.x += a * (meas.x - b_.x);
    b_.y += a * (meas.y - b_.y);
    b_.z += a * (meas.z - b_.z);
  }
private:
  float tau_ = 20.0f;
  Vec3f b_{0,0,0};
  bool inited_ = false;
};
