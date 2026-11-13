/* Built-in 2x12 cabinet and room impulse responses, synthesised from designed magnitude
 * curves. Port of makeCabIR / makeRoomIR in ../../index.html. Standard library only. */
#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

namespace bb {

inline void fft(std::vector<double>& re, std::vector<double>& im, bool inverse) {
  const size_t n = re.size();
  for (size_t i = 1, j = 0; i < n; i++) {
    size_t bit = n >> 1;
    for (; j & bit; bit >>= 1) j ^= bit;
    j ^= bit;
    if (i < j) { std::swap(re[i], re[j]); std::swap(im[i], im[j]); }
  }
  for (size_t len = 2; len <= n; len <<= 1) {
    const double ang = (inverse ? 2 : -2) * 3.141592653589793 / (double)len, wr = std::cos(ang), wi = std::sin(ang);
    for (size_t i = 0; i < n; i += len) {
      double cr = 1, ci = 0;
      for (size_t k = 0; k < len / 2; k++) {
        const size_t a = i + k, b = a + len / 2;
        const double xr = re[b] * cr - im[b] * ci, xi = re[b] * ci + im[b] * cr;
        re[b] = re[a] - xr; im[b] = im[a] - xi; re[a] += xr; im[a] += xi;
        const double t = cr * wr - ci * wi; ci = cr * wi + ci * wr; cr = t;
      }
    }
  }
  if (inverse) for (size_t i = 0; i < n; i++) { re[i] /= (double)n; im[i] /= (double)n; }
}

struct Rng {
  uint32_t s;
  explicit Rng(uint32_t seed) : s(seed) {}
  double operator()() { s = s * 1664525u + 1013904223u; return s / 4294967296.0; }
};

/* Minimum-phase impulse response from a magnitude curve in dB. */
template <class F>
inline std::vector<double> minPhaseIR(F magDb, size_t N, double fs) {
  std::vector<double> re(N, 0.0), im(N, 0.0);
  for (size_t k = 0; k <= N / 2; k++) {
    const double v = std::max(-90.0, magDb(std::max((double)k, 0.5) * fs / (double)N)) * 2.302585092994046 / 20;
    re[k] = v; if (k > 0 && k < N / 2) re[N - k] = v;
  }
  fft(re, im, true);
  for (size_t i = 1; i < N / 2; i++) { re[i] *= 2; re[N - i] = 0; }
  std::fill(im.begin(), im.end(), 0.0);
  fft(re, im, false);
  for (size_t i = 0; i < N; i++) { const double m = std::exp(re[i]), p = im[i]; re[i] = m * std::cos(p); im[i] = m * std::sin(p); }
  fft(re, im, true);
  return re;
}

/* 2x12 open-back cabinet with two slightly different 12" drivers (L / R), close mic'd.
 * pos: 0 = cap, 1 = cone edge. dist: 0 = on the grille, 1 = pulled back. Fills ir[0], ir[1]. */
inline void makeCabIR(double fs, double pos, double dist, bool ribbon, std::vector<float> ir[2]) {
  const size_t N = 8192, len = 4096;
  auto sq = [](double x) { return x * x; };
  auto hp = [&](double f, double fc, double q) { const double r = f / fc; return 20 * std::log10(r * r / std::sqrt(sq(1 - r * r) + sq(r / q))); };
  auto lp = [&](double f, double fc, double q) { const double r = f / fc; return -10 * std::log10(sq(1 - r * r) + sq(r / q)); };
  auto bell = [&](double f, double fc, double g, double q) { return g / (1 + sq(q * (f / fc - fc / f))); };
  auto step = [](double f, double a, double b) { const double t = std::min(1.0, std::max(0.0, std::log(f / a) / std::log(b / a))); return t * t * (3 - 2 * t); };
  for (int ch = 0; ch < 2; ch++) {
    Rng r(1962u + (uint32_t)ch * 77u);
    double modes[34][3];
    for (auto& m : modes) {
      m[0] = 420 * std::pow(22.0, r());
      const double a = r() * 2 - 1, b = 1 + 2.2 * r();
      m[1] = a * b;
      m[2] = 7 + 16 * r();
    }
    const double sk = ch ? 1.035 : 1, tau = ch ? 0.00215 : 0.0019, amt = 0.14 + 0.3 * dist;
    auto mag = [&](double f) {
      double d = hp(f, 96 * sk, 1.35) + bell(f, 185, 1.5, 1.1) + bell(f, 430, -3, 1.4) + bell(f, 900, -1.2, 1.5)
        + bell(f, 1500 * sk, 2.5, 1.2) + bell(f, 2500 * sk, 6 - 2.5 * pos, 1.6) + bell(f, 3500 * sk, 4, 3)
        + bell(f, 4550 * sk, -8, 6) + bell(f, 5300 * sk, 2.2, 6) + lp(f, 5600 * sk, 0.9) + lp(f, 7600, 0.7);
      for (const auto& m : modes) d += bell(f, m[0], m[1], m[2]);
      d += (2 - 9 * pos) * step(f, 1400, 6000);                       // cap -> cone edge
      d += (1 - dist) * 5 * (1 - step(f, 70, 260));                   // proximity effect
      const double a = amt / (1 + sq(f / 1200)), w = 2 * 3.141592653589793 * f * tau;  // rear wave / floor bounce
      d += 10 * std::log10(1 + a * a - 2 * a * std::cos(w));
      if (ribbon) d += bell(f, 150, 2, 0.8) - 4.5 * std::log2(std::max(1.0, f / 2500)) + lp(f, 15000, 0.7);
      else d += hp(f, 55, 0.7) + bell(f, 5600, 4, 1.4) + bell(f, 350, -1, 0.8) + lp(f, 13500, 0.7);
      return d;
    };
    double ref = 0; int cnt = 0;
    for (double f = 200; f < 2500; f *= 1.06) { ref += mag(f); cnt++; }
    ref /= cnt;
    const std::vector<double> h = minPhaseIR([&](double f) { return mag(f) - ref - 5; }, N, fs);
    ir[ch].resize(len);
    for (size_t i = 0; i < len; i++) ir[ch][i] = (float)(h[i] * (i > len - 600 ? (double)(len - i) / 600 : 1.0));
  }
}

/* Short, dark stereo room. Fills ir[0], ir[1]. */
inline void makeRoomIR(double fs, std::vector<float> ir[2]) {
  const size_t len = (size_t)std::lround(fs * 0.7);
  for (int ch = 0; ch < 2; ch++) {
    Rng r(40u + (uint32_t)ch);
    std::vector<float>& d = ir[ch];
    d.assign(len, 0.0f);
    const size_t pre = (size_t)std::lround(fs * (0.009 + ch * 0.002));
    double lpf = 0, hpf = 0, e = 0;
    for (size_t i = pre; i < len; i++) {
      const double t = (double)(i - pre) / fs, x = (r() * 2 - 1) * std::exp(-t / 0.085) * (1 - std::exp(-t / 0.004));
      lpf += 0.32 * (x - lpf); hpf += 0.02 * (lpf - hpf);
      d[i] = (float)(lpf - hpf); e += (double)d[i] * d[i];
    }
    const float g = (float)(0.6 / std::sqrt(e));
    for (size_t i = 0; i < len; i++) d[i] *= g;
  }
}

}  // namespace bb
