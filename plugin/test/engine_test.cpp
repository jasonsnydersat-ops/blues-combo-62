// Runs the amp engine on fixed test signals. With an argument, writes the output as raw
// float32 to that file (used to compare against the JavaScript engine); always checks sanity.
#include "../src/AmpCore.h"
#include <chrono>
#include <cstdio>

struct Scn { int mode; double sens; bool trem; double k[8]; double amp; double secs; };
static const Scn SCN[] = {
  { 0, 1, false, { 0.6, 0.0, 0.6, 0.6, 0.5, 0.5, 0.5, 0.0 }, 0.20, 0.5 },
  { 2, 1, false, { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.5, 0.0 }, 0.50, 0.5 },
  { 1, 1, true,  { 0.0, 0.4, 0.7, 0.5, 0.5, 0.5, 0.6, 0.8 }, 0.15, 1.5 },
};

int main(int argc, char** argv) {
  const double fs = 48000; const int os = 4;
  FILE* f = argc > 1 ? std::fopen(argv[1], "wb") : nullptr;
  bool ok = true;
  for (const Scn& s : SCN) {
    auto t0 = std::chrono::steady_clock::now();
    bb::AmpCore c(fs, os);
    auto t1 = std::chrono::steady_clock::now();
    c.mode = s.mode; c.sens = s.sens; c.setTremolo(s.trem);
    for (int k = 0; k < 8; k++) c.setKnob((bb::AmpCore::Knob)k, s.k[k]);
    const int n = (int)(fs * s.secs);
    std::vector<float> in((size_t)n), out((size_t)n);
    for (int i = 0; i < n; i++) {
      const double t = i / fs;
      in[(size_t)i] = (float)(s.amp * std::exp(-t * 3) * (std::sin(2 * 3.141592653589793 * 110 * t) + 0.5 * std::sin(2 * 3.141592653589793 * 329.6 * t)));
    }
    for (int i = 0; i < n; i += 128) c.process(&in[(size_t)i], &out[(size_t)i], std::min(128, n - i));
    auto t2 = std::chrono::steady_clock::now();
    double pk = 0, e = 0;
    for (float v : out) { if (!(v == v)) ok = false; pk = std::max(pk, (double)std::fabs(v)); e += v * v; }
    if (!(pk > 1e-3 && pk < 4)) ok = false;
    std::printf("build %.0f ms, run %.0f ms for %.1f s audio (%.1f%% of real time), peak %.4f, rms %.4f, HT %.1f V\n",
      std::chrono::duration<double, std::milli>(t1 - t0).count(), std::chrono::duration<double, std::milli>(t2 - t1).count(), s.secs,
      100 * std::chrono::duration<double>(t2 - t1).count() / s.secs, pk, std::sqrt(e / n), c.Vb);
    if (f) std::fwrite(out.data(), sizeof(float), out.size(), f);
  }
  if (f) std::fclose(f);
  std::puts(ok ? "OK" : "FAILED");
  return ok ? 0 : 1;
}
