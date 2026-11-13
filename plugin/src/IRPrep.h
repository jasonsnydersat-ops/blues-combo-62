/* Prepares a user-supplied cabinet impulse response: resample to the host rate, strip
 * leading/trailing silence, cap the length, fade the tail and level-match it to the
 * built-in cabinet. Mirrors prepIR() in ../../index.html. Standard library only. */
#pragma once

#include <algorithm>
#include <cmath>
#include <vector>

namespace bb {

/* Band-limited (Kaiser-windowed sinc) resampling of one channel. */
inline std::vector<float> resampleIR(const std::vector<float>& in, double fsIn, double fsOut) {
  if (in.empty() || fsIn == fsOut) return in;
  const double ratio = fsOut / fsIn, cut = std::min(1.0, ratio) * 0.96, pi = 3.141592653589793;
  const int half = (int)std::ceil(24 / cut);
  const double beta = 9.0;
  auto i0 = [](double x) { double s = 1, t = 1; for (int k = 1; k < 30; k++) { t *= (x / (2 * k)) * (x / (2 * k)); s += t; } return s; };
  const double i0b = i0(beta);
  const size_t nOut = (size_t)std::floor((double)in.size() * ratio);
  std::vector<float> out(nOut);
  for (size_t i = 0; i < nOut; i++) {
    const double t = (double)i / ratio;
    const long c = (long)std::floor(t);
    double acc = 0;
    for (long k = c - half + 1; k <= c + half; k++) {
      if (k < 0 || k >= (long)in.size()) continue;
      const double d = t - (double)k, r = d / half;
      if (r <= -1 || r >= 1) continue;
      const double w = i0(beta * std::sqrt(1 - r * r)) / i0b;
      const double x = pi * d * cut;
      acc += in[(size_t)k] * (std::fabs(x) < 1e-9 ? 1.0 : std::sin(x) / x) * cut * w;
    }
    out[i] = (float)acc;
  }
  return out;
}

/* ch: 1 or 2 channels at sample rate fs (already at the host rate). Returns false if silent.
 * level is the target root-sum-square of the louder channel. */
inline bool trimAndLevelIR(std::vector<std::vector<float>>& ch, double fs, double level = 0.45, double maxSeconds = 1.5) {
  if (ch.empty() || ch[0].empty()) return false;
  if (ch.size() > 2) ch.resize(2);
  const size_t len = ch[0].size();
  float pk = 0;
  for (const auto& d : ch) for (float v : d) pk = std::max(pk, std::fabs(v));
  if (!(pk > 0)) return false;
  size_t a = len, b = 0;
  for (const auto& d : ch) {
    size_t i = 0; while (i < d.size() && std::fabs(d[i]) < pk * 1e-3f) i++;
    size_t j = d.size() - 1; while (j > i && std::fabs(d[j]) < pk * 1e-4f) j--;
    a = std::min(a, i); b = std::max(b, j);
  }
  a = a > 4 ? a - 4 : 0;
  const size_t n = std::max<size_t>(64, std::min<size_t>(b - a + 1, (size_t)std::lround(fs * maxSeconds)));
  const size_t fade = std::min(n / 2, (size_t)std::lround(fs * 0.004));
  double en = 0;
  for (auto& d : ch) {
    std::vector<float> o(n);
    double e = 0;
    for (size_t i = 0; i < n; i++) {
      const float src = a + i < d.size() ? d[a + i] : 0.0f;
      const float x = src * (i >= n - fade ? (float)(n - 1 - i) / (float)fade : 1.0f);
      o[i] = x; e += (double)x * x;
    }
    d.swap(o);
    en = std::max(en, e);
  }
  const float g = (float)(level / std::sqrt(en));
  for (auto& d : ch) for (float& v : d) v *= g;
  return true;
}

}  // namespace bb
