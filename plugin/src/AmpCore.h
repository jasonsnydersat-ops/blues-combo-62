/* 1962-style 2x12 combo amplifier model: circuit-level DSP core.
 *
 * This is a line-for-line C++ port of ../../amp-processor.js (the web version's engine).
 * It has no dependencies beyond the C++ standard library.
 *
 * Signal path (all oversampled):
 *   input jacks -> V1 (ECC83, two halves, shared cathode) -> loudness pots + mixers
 *   -> V2a gain stage DC-coupled to V2b cathode follower -> tone stack
 *   -> long-tail phase inverter -> 0.1uF couplers with grid conduction
 *   -> push-pull KT66 into a reactive speaker load through the output transformer,
 *      negative feedback + presence, solved implicitly every sample
 *   -> GZ34 rectifier + reservoir + choke-fed screens (sag and ripple)
 *   Tremolo: V6a phase-shift oscillator -> V6b follower -> J174 FET shunting channel II's
 *   mixing network.
 *
 * Nothing in process() or the knob updates allocates memory.
 */
#pragma once

#include <algorithm>
#include <cmath>
#include <cstring>
#include <vector>

namespace bb {

/* ------------------------------------------------------------------ tables */
struct LUT {
  std::vector<double> d;
  double x0 = 0, x1 = 0, inv = 0;
  int n = 0;
};

template <class F>
inline void fillLUT(LUT& t, int n, double x0, double x1, F fn) {
  t.d.resize((size_t)n + 1);
  const double dx = (x1 - x0) / (n - 1);
  for (int i = 0; i < n; i++) t.d[(size_t)i] = fn(x0 + i * dx);
  t.d[(size_t)n] = t.d[(size_t)n - 1];
  t.x0 = x0; t.x1 = x1; t.inv = 1 / dx; t.n = n;
}
inline double lut(const LUT& t, double x) {
  const double p = (x - t.x0) * t.inv;
  if (p <= 0) return t.d[0];
  const int m = t.n - 1;
  if (p >= m) return t.d[(size_t)m];
  const int i = (int)p;
  const double a = t.d[(size_t)i];
  return a + (t.d[(size_t)i + 1] - a) * (p - i);
}
/* as lut(), also returning the local slope */
inline double lutS(const LUT& t, double x, double& sl) {
  const double p = (x - t.x0) * t.inv;
  if (p <= 0) { sl = 0; return t.d[0]; }
  const int m = t.n - 1;
  if (p >= m) { sl = 0; return t.d[(size_t)m]; }
  const int i = (int)p;
  const double a = t.d[(size_t)i], s = t.d[(size_t)i + 1] - a;
  sl = s * t.inv;
  return a + s * (p - i);
}

/* Inverse of m(v) = v + Rp*ig(v): the grid voltage reached when an open-circuit
 * voltage x drives a conducting grid through resistance Rp. Identity below vA.
 * Rebuilding into an existing table of the same size does not allocate. */
template <class F>
inline void buildClamp(LUT& t, F ig, double Rp, double lo, double hi, int n, double vA, double vB) {
  enum { M = 2048 };
  double mv[M + 1];
  const double dv = (vB - vA) / M;
  for (int j = 0; j <= M; j++) { const double v = vA + j * dv; mv[j] = v + Rp * ig(v); }
  if ((int)t.d.size() != n + 1) t.d.resize((size_t)n + 1);
  t.x0 = lo; t.x1 = hi; t.inv = (n - 1) / (hi - lo); t.n = n;
  const double dx = (hi - lo) / (n - 1);
  double* d = t.d.data();
  const double m0 = mv[0], endSlope = dv / (mv[M] - mv[M - 1]);
  int j = 0;
  for (int i = 0; i < n; i++) {
    const double x = lo + i * dx;
    if (x <= m0) { d[i] = x - (m0 - vA); continue; }
    while (j < M - 1 && mv[j + 1] < x) j++;
    if (x >= mv[M]) { d[i] = vB + (x - mv[M]) * endSlope; continue; }
    d[i] = vA + (j + (x - mv[j]) / (mv[j + 1] - mv[j])) * dv;
  }
  d[n] = d[n - 1];
}

/* ------------------------------------------------------------ valve models */
struct Valve { double mu, ex, kg1, kg2, kp, kvb; };
static const Valve ECC83 = { 100, 1.4, 1060, 0, 600, 300 };
static const Valve KT66 = { 8.2, 1.35, 1350, 4500, 45, 13 };

inline double triIp(double vgk, double vpk, const Valve& T) {
  if (vpk <= 0) return 0;
  const double t = T.kp * (1 / T.mu + vgk / std::sqrt(T.kvb + vpk * vpk));
  const double e1 = (vpk / T.kp) * (t > 30 ? t : std::log1p(std::exp(t)));
  return 2 * std::pow(e1, T.ex) / T.kg1;
}
inline double gridI(double vgk, double g, double s, double vo) {
  const double x = (vgk + vo) / s;
  const double sp = x > 30 ? vgk + vo : s * std::log1p(std::exp(x));
  return g * sp * std::sqrt(sp);
}
inline double igT(double v) { return gridI(v, 5e-4, 0.12, 0.2); }   // ECC83 grid
inline double igP(double v) { return gridI(v, 1.5e-3, 0.3, 0.0); }  // KT66 grid

template <class F>
inline double bisect(F f, double lo, double hi, int n) {  // f decreasing; returns root
  for (int i = 0; i < n; i++) { const double m = 0.5 * (lo + hi); if (f(m) > 0) lo = m; else hi = m; }
  return 0.5 * (lo + hi);
}
/* plate current of a triode with plate load Rp from supply Vs (cathode fixed) */
inline double loadLine(double vgk, double Vs, double Rp, const Valve& T) {
  if (Vs <= 0) return 0;
  return bisect([&](double I) { return triIp(vgk, Vs - I * Rp, T) - I; }, 0, Vs / Rp, 26);
}

/* --------------------------------------------- netlist -> discrete state space */
struct El { char t; int a, b; double v; };   // 'R' | 'C' | 'L'; node 0 is ground
struct In { char t; int a; };                // 'V' | 'I'

enum { SS_MAXS = 5, SS_MAXI = 3, SS_MAXN = 20, SS_MAXE = SS_MAXS + SS_MAXI };

struct SS {
  int ns = 0, ni = 0;
  double P[SS_MAXS * SS_MAXS] = {}, Q[SS_MAXS * SS_MAXI] = {}, C[SS_MAXS] = {}, D[SS_MAXI] = {};
  double x[SS_MAXS] = {}, z[SS_MAXS] = {}, up[SS_MAXI] = {};
  void clearState() { for (double& v : x) v = 0; for (double& v : up) v = 0; }
};

/* Gauss-Jordan with partial pivoting; A is N x N, B is N x ne. The solution replaces B. */
template <int RA, int CA, int CB>
inline void solveLin(double A[RA][CA], double B[RA][CB], int N, int ne) {
  for (int c = 0; c < N; c++) {
    int p = c; double best = std::fabs(A[c][c]);
    for (int r = c + 1; r < N; r++) { const double v = std::fabs(A[r][c]); if (v > best) { best = v; p = r; } }
    if (p != c) {
      for (int k = 0; k < N; k++) std::swap(A[p][k], A[c][k]);
      for (int k = 0; k < ne; k++) std::swap(B[p][k], B[c][k]);
    }
    const double piv = 1 / A[c][c];
    for (int k = c; k < N; k++) A[c][k] *= piv;
    for (int k = 0; k < ne; k++) B[c][k] *= piv;
    for (int r = 0; r < N; r++) {
      if (r == c) continue;
      const double f = A[r][c];
      if (f == 0) continue;
      for (int k = c; k < N; k++) A[r][k] -= f * A[c][k];
      for (int k = 0; k < ne; k++) B[r][k] -= f * B[c][k];
    }
  }
}

/* States: capacitor voltages (a-b) then inductor currents (a->b). Keeps ss.x / ss.up. */
inline void buildSS(SS& ss, int nNodes, const El* els, int nEls, const In* ins, int ni, int outNode, double fs) {
  const El* caps[SS_MAXS]; const El* inds[SS_MAXS];
  int nc = 0, nl = 0;
  for (int i = 0; i < nEls; i++) { if (els[i].t == 'C') caps[nc++] = &els[i]; else if (els[i].t == 'L') inds[nl++] = &els[i]; }
  const int ns = nc + nl, ne = ns + ni;
  int nvs = nc; for (int i = 0; i < ni; i++) if (ins[i].t == 'V') nvs++;
  const int N = nNodes + nvs;
  double A[SS_MAXN][SS_MAXN], B[SS_MAXN][SS_MAXE];
  for (int i = 0; i < N; i++) { for (int j = 0; j < N; j++) A[i][j] = 0; for (int j = 0; j < ne; j++) B[i][j] = 0; }
  for (int i = 0; i < nEls; i++) {
    const El& e = els[i];
    if (e.t != 'R') continue;
    const double g = 1 / std::max(e.v, 1.0); const int a = e.a - 1, b = e.b - 1;
    if (a >= 0) A[a][a] += g;
    if (b >= 0) A[b][b] += g;
    if (a >= 0 && b >= 0) { A[a][b] -= g; A[b][a] -= g; }
  }
  int row = nNodes;
  auto vsrc = [&](int a, int b, int exc) {
    if (a >= 0) { A[row][a] = 1; A[a][row] = 1; }
    if (b >= 0) { A[row][b] = -1; A[b][row] = -1; }
    B[row][exc] = 1; row++;
  };
  for (int k = 0; k < nc; k++) vsrc(caps[k]->a - 1, caps[k]->b - 1, k);
  for (int k = 0; k < nl; k++) {
    if (inds[k]->a > 0) B[inds[k]->a - 1][nc + k] -= 1;
    if (inds[k]->b > 0) B[inds[k]->b - 1][nc + k] += 1;
  }
  for (int k = 0; k < ni; k++) {
    if (ins[k].t == 'V') vsrc(ins[k].a - 1, -1, ns + k); else B[ins[k].a - 1][ns + k] += 1;
  }
  solveLin<SS_MAXN, SS_MAXN, SS_MAXE>(A, B, N, ne);   // B is now X
  double Mx[SS_MAXS][SS_MAXE];
  for (int k = 0; k < nc; k++) for (int j = 0; j < ne; j++) Mx[k][j] = B[nNodes + k][j] / caps[k]->v;
  for (int k = 0; k < nl; k++) {
    const El& e = *inds[k];
    for (int j = 0; j < ne; j++) Mx[nc + k][j] = ((e.a > 0 ? B[e.a - 1][j] : 0) - (e.b > 0 ? B[e.b - 1][j] : 0)) / e.v;
  }
  const double h2 = 0.5 / fs;
  double L[SS_MAXS][SS_MAXS], R[SS_MAXS][SS_MAXE];
  for (int i = 0; i < ns; i++) {
    for (int j = 0; j < ns; j++) { L[i][j] = (i == j ? 1 : 0) - h2 * Mx[i][j]; R[i][j] = (i == j ? 1 : 0) + h2 * Mx[i][j]; }
    for (int j = 0; j < ni; j++) R[i][ns + j] = h2 * Mx[i][ns + j];
  }
  solveLin<SS_MAXS, SS_MAXS, SS_MAXE>(L, R, ns, ns + ni);   // R is now S
  ss.ns = ns; ss.ni = ni;
  for (int i = 0; i < ns; i++) {
    for (int j = 0; j < ns; j++) ss.P[i * ns + j] = R[i][j];
    for (int j = 0; j < ni; j++) ss.Q[i * ni + j] = R[i][ns + j];
    ss.C[i] = B[outNode - 1][i];
  }
  for (int j = 0; j < ni; j++) ss.D[j] = B[outNode - 1][ns + j];
}
/* single-input step: x <- P x + Q (u_prev + u), y = C x + D u */
inline double ssStep(SS& s, double u) {
  const int ns = s.ns; const double uu = s.up[0] + u;
  double y = s.D[0] * u;
  for (int i = 0; i < ns; i++) {
    double a = s.Q[i] * uu; const int o = i * ns;
    for (int j = 0; j < ns; j++) a += s.P[o + j] * s.x[j];
    s.z[i] = a; y += s.C[i] * a;
  }
  for (int i = 0; i < ns; i++) s.x[i] = s.z[i];
  s.up[0] = u;
  return y;
}

/* ------------------------------------------------------------- oversampling */
inline std::vector<double> kaiserSinc(int taps, double fc, double beta) {
  std::vector<double> h((size_t)taps);
  const double mid = (taps - 1) / 2.0, pi = 3.141592653589793;
  auto i0 = [](double x) { double s = 1, t = 1; for (int k = 1; k < 30; k++) { t *= (x / (2 * k)) * (x / (2 * k)); s += t; } return s; };
  double sum = 0;
  for (int i = 0; i < taps; i++) {
    const double n = i - mid, r = n / mid;
    const double w = i0(beta * std::sqrt(std::max(0.0, 1 - r * r))) / i0(beta);
    h[(size_t)i] = (n == 0 ? 2 * fc : std::sin(2 * pi * fc * n) / (pi * n)) * w;
    sum += h[(size_t)i];
  }
  for (int i = 0; i < taps; i++) h[(size_t)i] /= sum;
  return h;
}

/* Tremolo, values from the 1962 reissue schematic (V6 ECC83 + J174). */
namespace OSC { const double Vs = 320, RL = 150e3, Rk = 560, Rfs = 2200, Cpk = 47e-9, C1 = 32e-9, R1 = 68e3, Rpot = 2e6, C2 = 22e-9, R2 = 503e3, C3 = 22e-9, R3 = 1.47e6; }
namespace FET { const double vp = 7, idss = 0.045; }   // J174, P-channel, source grounded
inline double fetI(double v, double vg) {              // current from the drain node into the FET
  const double beta = FET::idss / (FET::vp * FET::vp);
  if (v <= 0) {
    const double vov = FET::vp - vg, vds = -v;
    if (vov <= 0) return 0;
    return -(vds < vov ? beta * (2 * vov - vds) * vds : beta * vov * vov);
  }
  const double vov = FET::vp - vg + v;
  if (vov <= 0) return 0;
  return v < vov ? beta * (2 * vov - v) * v : beta * vov * vov;
}

inline double audioTaper(double x) { return (std::pow(81.0, x) - 1) / 80; }  // 10% at half rotation

/* ================================================================ amp core */
class AmpCore {
public:
  enum Knob { VOL1, VOL2, TREBLE, MID, BASS, PRESENCE, SPEED, DEPTH, NUM_KNOBS };

  int mode = 0;          // 0 high treble (channel I), 1 normal (channel II), 2 jumpered
  double sens = 1;       // 1 = high-sensitivity jack, 0.5 = low
  double inScale = 2.0;  // volts at digital full scale
  double Vb = 445;       // HT on the plates (for a meter)
  long iters = 0;

  AmpCore(double sampleRate, int oversampling)
    : fs(sampleRate), os(oversampling), fo(sampleRate * oversampling), dt(1 / (sampleRate * oversampling)) {
    h = kaiserSinc(per * os + 1, 0.47 / os, 7.9);
    hLen = (int)h.size();

    buildValves();

    const double p0[NUM_KNOBS] = { 0.6, 0.0, 0.6, 0.6, 0.5, 0.5, 0.5, 0 };
    for (int i = 0; i < NUM_KNOBS; i++) p[i] = tgt[i] = p0[i];

    vk1 = vk1Q;
    sc = std::pow(Vg2 / KT66.kp, KT66.ex) * 2 / KT66.kg1;
    hpC = 1 - std::exp(-2 * 3.141592653589793 * 7 / fo);

    buildLoad();
    updateNetworks();
    calibrateSupply();
    // Let the tremolo oscillator build up at its slowest setting (smallest swing) and set the gate
    // trimmer (RV2) so the LFO peak still pinches the FET off there; then return to the knob setting.
    const double spd = p[SPEED]; p[SPEED] = 0; updateNetworks();
    osc.x[0] = 60;
    const int nfs = (int)fs;
    for (int i = 0; i < nfs * 5; i++) stepLFO();
    double top = 0; for (int i = 0; i < nfs; i++) { stepLFO(); if (oscVp > top) top = oscVp; }
    fetK = (FET::vp + 0.6) / top;
    p[SPEED] = spd; updateNetworks();
    for (int i = 0; i < nfs * 2; i++) stepLFO();
    float z[256], o[256];
    const int warm = (int)std::ceil(fs * 0.3 / 256);
    for (int i = 0; i < warm; i++) { std::memset(z, 0, sizeof z); process(z, o, 256); }
  }

  /* knob values are 0..1 */
  void setKnob(Knob k, double v) { tgt[k] = std::min(1.0, std::max(0.0, v)); }
  double knob(Knob k) const { return tgt[k]; }
  void setTremolo(bool on) { if (on != tremOn) { tremOn = on; dirty = true; } }
  /* latency of the oversampling filters, in base-rate samples */
  int latencySamples() const { return per; }

  /* inp/out: n samples at the base rate (may be the same buffer). Output is normalised half-primary voltage. */
  void process(const float* inp, float* out, int n) {
    smooth();
    const Valve& PV = KT66;
    const double ikvb = 1 / PV.kvb, kgr = PV.kg1 / PV.kg2;
    const double kCk = dt / Ck1, iRk = 1 / Rk1;
    const double gB = mode == 1 ? 0 : sens * inScale, gN = mode == 0 ? 0 : sens * inScale;
    const double piRsI = 1 / piRs, piCk = dt / piC;
    const double pwRsI = 1 / pwRs, pwCk = dt / pwC;
    double* mx = mix.x; double* mz = mix.z; const double* mP = mix.P; const double* mQ = mix.Q; const double* mC = mix.C; const double* mD = mix.D;
    double* mup = mix.up; const double* dC = mixD.C; const double* dD = mixD.D;
    const double Zm = fetZm, Zd = fetZd;
    const double acInc = 2 * 3.141592653589793 * 50 * dt, iRr = 1 / Rrect, iRch = 1 / Rch, kC1 = dt / C1, kC2 = dt / C2;
    double* lx = load.x; double* lz = load.z; const double* lP = load.P; const double* lQ = load.Q; const double* lC = load.C; const int lns = load.ns;
    const double fP = fb.P[0], fQ = fb.Q[0], fC = fb.C[0];
    const double* hh = h.data();
    long it_count = 0;

    for (int s = 0; s < n; s++) {
      upPos = (upPos + 1) & 63; upBuf[upPos] = inp[s];
      stepLFO();
      fetG = fetK * (oscVp - 1);   // V6b follower -> R40 / RV2 divider -> J174 gate
      for (int ph = 0; ph < os; ph++) {
        // ---- polyphase interpolation
        double x = 0;
        for (int k = 0, j = ph; k < per; k++, j += os) x += hh[j] * upBuf[(upPos - k) & 63];
        x *= os;

        // ---- V1
        const double ia = lut(tabV1, lut(clampIn, x * gB - vk1));
        const double ib = lut(tabV1, lut(clampIn, x * gN - vk1));
        vk1 += kCk * (ia + ib - vk1 * iRk);
        double u1 = -Rp1 * (ia - ip1Q), u2 = -Rp1 * (ib - ip1Q);
        hp1 += hpC * (u1 - hp1); u1 -= hp1;
        hp2 += hpC * (u2 - hp2); u2 -= hp2;

        // ---- loudness pots + mixers, V2a + follower, tone stack
        // (the tremolo FET is solved implicitly against the node it shunts)
        const double b0 = mup[0] + u1, b1 = mup[1] + u2, b2 = mup[2];
        double ym = mD[0] * u1 + mD[1] * u2, yd = dD[0] * u1 + dD[1] * u2;
        for (int i = 0; i < 5; i++) {
          const int k = 5 * i, q = 3 * i;
          const double a = mQ[q] * b0 + mQ[q + 1] * b1 + mQ[q + 2] * b2 + mP[k] * mx[0] + mP[k + 1] * mx[1] + mP[k + 2] * mx[2] + mP[k + 3] * mx[3] + mP[k + 4] * mx[4];
          mz[i] = a; ym += mC[i] * a; yd += dC[i] * a;
        }
        for (int it = 0; it < 6; it++) {
          const double i0 = fetI(fetV, fetG), di = (fetI(fetV + 1e-3, fetG) - i0) * 1e3;
          const double stp = (fetV - yd + Zd * i0) / (1 + Zd * di);
          fetV -= stp;
          if (stp < 1e-6 && stp > -1e-6) break;
        }
        const double u3 = (fetV - yd) / Zd;
        for (int i = 0; i < 5; i++) mx[i] = mz[i] + mQ[3 * i + 2] * u3;
        mup[0] = u1; mup[1] = u2; mup[2] = u3;
        const double vm = ym + Zm * u3;
        const double v2 = lut(tabV2, lut(clampV2, vm));
        const double vt = ssStep(tone, v2);

        // ---- phase inverter input coupler (grid conduction charges the 22n)
        const double vn = lut(clampPI, piK * (vt - cPI) - vF) + vF;
        cPI += piCk * (vt - cPI - vn) * piRsI;

        // ---- power stage: implicit solve for the half-primary voltage
        double Veq = 0;
        for (int i = 0; i < lns; i++) {
          double a = lQ[i] * uL; const int o = i * lns;
          for (int j = 0; j < lns; j++) a += lP[o + j] * lx[j];
          lz[i] = a; Veq += lC[i] * a;
        }
        const double zF = fP * fb.x[0] + fQ * nSec * VhPrev, hF = fC * zF;
        double lo = -1e5, hi = 1e5, va = 0, vb = 0, gA = 0, gB2 = 0, IA = 0, IB = 0, LA = 0, LB = 0, SL = 0;
        const double iVg2 = 1 / Vg2;
        for (int it = 0; it < 14; it++) {
          it_count++;
          vF = aF * Vh + hF;
          const double vd = vn - vF;
          va = lutS(tabPIa, vd, SL); const double sa = SL;
          vb = lutS(tabPIb, vd, SL); const double sb = SL;
          gA = lutS(clampPw, pwK * (vb - cA), SL); const double sgA = SL * pwK * sb;
          gB2 = lutS(clampPw, pwK * (va - cB), SL); const double sgB = SL * pwK * sa;
          LA = lutS(tabL, (vBias + gA) * iVg2, SL); const double sLA = SL * iVg2 * sgA;
          LB = lutS(tabL, (vBias + gB2) * iVg2, SL); const double sLB = SL * iVg2 * sgB;
          const double pA = (Vb - Vh) * ikvb, pB = (Vb + Vh) * ikvb;
          const double atA = pA > 0 ? std::atan(pA) : 0, atB = pB > 0 ? std::atan(pB) : 0;
          const double dA = pA > 0 ? ikvb / (1 + pA * pA) : 0, dB = pB > 0 ? ikvb / (1 + pB * pB) : 0;
          IA = sc * LA * atA; IB = sc * LB * atB;
          const double g = Vh - Req * (IA - IB) - Veq;
          const double dg = 1 - Req * sc * ((-sLA * aF * atA - LA * dA) - (-sLB * aF * atB + LB * dB));
          if (g > 0) hi = Vh; else lo = Vh;
          double nx = Vh - g / dg;
          if (!(nx > lo && nx < hi)) nx = 0.5 * (lo + hi);
          const double e = nx - Vh; Vh = nx;
          if (e < 2e-3 && e > -2e-3) break;
        }
        vF = aF * Vh + hF;
        const double u = (Vh - Veq) / Req;
        for (int i = 0; i < lns; i++) lx[i] = lz[i] + lQ[i] * u;
        uL = u;
        fb.x[0] = zF + fQ * nSec * Vh; VhPrev = Vh;
        cA += pwCk * (vb - cA - gA) * pwRsI;
        cB += pwCk * (va - cB - gB2) * pwRsI;

        // ---- GZ34 + reservoir, choke-fed screens
        acPh += acInc; if (acPh > 6.283185307) acPh -= 6.283185307;
        double vac = Vpk * std::sin(acPh); if (vac < 0) vac = -vac;
        const double ich = (Vb - Vg2) * iRch;
        Vb += kC1 * ((vac > Vb ? (vac - Vb) * iRr : 0) - IA - IB - ich);
        Vg2 += kC2 * (ich - sc * (LA + LB) * kgr - iPre);
        sc = std::pow(Vg2 / PV.kp, PV.ex) * 2 / PV.kg1;

        dnPos = (dnPos + 1) & 1023; dnBuf[dnPos] = Vh;
      }
      // ---- decimation
      double y = 0;
      for (int j = 0; j < hLen; j++) y += hh[j] * dnBuf[(dnPos - j) & 1023];
      out[s] = (float)(y * (1.0 / 430));
    }
    iters = it_count;
    if (!(Vh == Vh) || !(vk1 == vk1) || !(fetV == fetV) || !(Vb == Vb)) resetState();
  }

  void resetState() {
    vk1 = vk1Q; hp1 = hp2 = 0; cPI = cA = cB = 0; vF = 0;
    Vh = 0; uL = 0; VhPrev = 0; Vb = 445; Vg2 = 438;
    fetV = 0;
    sc = std::pow(Vg2 / KT66.kp, KT66.ex) * 2 / KT66.kg1;
    mix.clearState(); tone.clearState(); fb.clearState(); load.clearState();
    std::memset(upBuf, 0, sizeof upBuf); std::memset(dnBuf, 0, sizeof dnBuf);
  }

private:
  const double fs; const int os; const double fo, dt;
  enum { per = 20 };
  std::vector<double> h; int hLen = 0;
  double upBuf[64] = {}; int upPos = 0;
  double dnBuf[1024] = {}; int dnPos = 0;

  double p[NUM_KNOBS], tgt[NUM_KNOBS];
  bool dirty = true, tremOn = true;

  // valve tables
  LUT tabV1, clampIn, tabV2, tabIg2, clampV2, tabPIa, tabPIb, clampPI, tabL, clampPw, clampOsc;
  bool haveClampV2 = false;
  double Rp1 = 0, Rk1 = 0, Ck1 = 0, vk1Q = 0, ip1Q = 0;
  double piRs = 0, piRl = 0, piC = 0, piK = 0;
  double vBias = 0, idle = 0, pwRs = 0, pwRl = 0, pwC = 0, pwK = 0;
  double rth = 0;

  // linear networks
  SS load, mix, mixD, tone, fb, osc, oscK, oscP;
  double Req = 0, nSec = 0, aF = 0, fetZm = 0, fetZd = 0, oscRz = 0;

  // supply
  double Rrect = 90, C1 = 32e-6, Rch = 250, C2 = 32e-6, iPre = 0.006, Vpk = 480;

  // dynamic state
  double vk1 = 0, hp1 = 0, hp2 = 0, cPI = 0, cA = 0, cB = 0, vF = 0;
  double Vh = 0, uL = 0, VhPrev = 0, Vg2 = 438, acPh = 0, sc = 0, hpC = 0;
  double oscIp = 0.0005, oscVp = 200, fetK = 0.035, fetV = 0, fetG = 20;

  void buildValves() {
    const Valve& T = ECC83;
    /* ---- V1: both halves, 100k plates, shared 820R / 250uF cathode ---- */
    const double Vb1 = 265, Rp = 100e3, Rk = 820;
    Rp1 = Rp; Rk1 = Rk; Ck1 = 250e-6;
    vk1Q = bisect([&](double vk) { return 2 * loadLine(-vk, Vb1 - vk, Rp, T) - vk / Rk; }, 0, 6, 30);
    fillLUT(tabV1, 4096, -9, 4, [&](double v) { return loadLine(v, Vb1 - vk1Q, Rp, T); });
    ip1Q = lut(tabV1, -vk1Q);
    buildClamp(clampIn, igT, 68e3, -20, 20, 8192, -3, 4);

    /* ---- V2a (100k / 820R unbypassed) DC-coupled to V2b follower (100k) ---- */
    const double Vb2 = 300, Rk2 = 820, Rcf = 100e3;
    auto solve = [&](double vg, double& outVk2, double& outIg) {
      double i1 = 0, vk2 = 0;
      bisect([&](double Vp) {
        i1 = bisect([&](double I) { return triIp(vg - Rk2 * I, Vp - Rk2 * I, T) - I; }, 0, 0.012, 24);
        vk2 = bisect([&](double Vk) { return triIp(Vp - Vk, Vb2 - Vk, T) + igT(Vp - Vk) - Vk / Rcf; }, 0, Vb2, 24);
        return (Vb2 - Vp) / Rp - i1 - igT(Vp - vk2);
      }, 0, Vb2, 24);
      outVk2 = vk2; outIg = igT(vg - Rk2 * i1);
    };
    const int n2 = 2048;
    double q2, dummy; solve(0, q2, dummy);
    fillLUT(tabV2, n2, -10, 8, [](double) { return 0.0; });
    fillLUT(tabIg2, n2, -10, 8, [](double) { return 0.0; });
    for (int i = 0; i < n2; i++) {
      double c, g; solve(-10 + i * 18.0 / (n2 - 1), c, g);
      tabV2.d[(size_t)i] = c - q2; tabIg2.d[(size_t)i] = g;
    }
    tabV2.d[n2] = tabV2.d[n2 - 1]; tabIg2.d[n2] = tabIg2.d[n2 - 1];

    /* ---- long-tail phase inverter: 82k / 100k, 470R bias, 10k + 5k tail ---- */
    const double Vb3 = 320, Ra = 82e3, Rb = 100e3, Rbias = 470, Rtail = 15e3, Rkk = Rbias + Rtail;
    const double itQ = bisect([&](double It) {
      return loadLine(-Rbias * It, Vb3 - Rkk * It, Ra, T) + loadLine(-Rbias * It, Vb3 - Rkk * It, Rb, T) - It;
    }, 0, 0.006, 30);
    const double vxQ = itQ * Rtail;
    auto piSolve = [&](double vd, double r[3]) {
      double ia = 0, ib = 0, ig = 0;
      bisect([&](double Vk) {
        ia = loadLine(vxQ + vd - Vk, Vb3 - Vk, Ra, T);
        ib = loadLine(vxQ - Vk, Vb3 - Vk, Rb, T);
        ig = igT(vxQ + vd - Vk);
        return ia + ib + ig - Vk / Rkk;
      }, 0, Vb3, 26);
      r[0] = ia; r[1] = ib; r[2] = ig;
    };
    const int n3 = 1024; double q3[3]; piSolve(0, q3);
    LUT tig;
    fillLUT(tabPIa, n3, -30, 30, [](double) { return 0.0; });
    fillLUT(tabPIb, n3, -30, 30, [](double) { return 0.0; });
    fillLUT(tig, n3, -30, 30, [](double) { return 0.0; });
    for (int i = 0; i < n3; i++) {
      double r[3]; piSolve(-30 + i * 60.0 / (n3 - 1), r);
      tabPIa.d[(size_t)i] = -Ra * (r[0] - q3[0]); tabPIb.d[(size_t)i] = -Rb * (r[1] - q3[1]); tig.d[(size_t)i] = r[2];
    }
    tabPIa.d[n3] = tabPIa.d[n3 - 1]; tabPIb.d[n3] = tabPIb.d[n3 - 1]; tig.d[n3] = tig.d[n3 - 1];
    // tone stack -> 22n -> PI grid, 1M leak
    piRs = 120e3; piRl = 1e6; piC = 22e-9;
    piK = piRl / (piRs + piRl);
    buildClamp(clampPI, [&](double v) { return lut(tig, v); }, piRs * piRl / (piRs + piRl), -150, 150, 8192, -30, 30);

    /* ---- KT66 pair ---- */
    const Valve& P = KT66;
    fillLUT(tabL, 4096, -0.5, 0.08, [&](double x) {
      const double t = P.kp * (1 / P.mu + x);
      return std::pow(t > 30 ? t : std::log1p(std::exp(t)), P.ex);
    });
    const double idleI = 0.036, VbN = 445, Vg2N = 438;
    const double scN = std::pow(Vg2N / P.kp, P.ex) * 2 / P.kg1;
    vBias = bisect([&](double vb) { return idleI - scN * lut(tabL, vb / Vg2N) * std::atan(VbN / P.kvb); }, -90, -10, 40);
    idle = idleI;
    pwRs = 50e3; pwRl = 220e3; pwC = 0.1e-6;
    pwK = pwRl / (pwRs + pwRl);
    const double vb = vBias;
    buildClamp(clampPw, [vb](double v) { return igP(vb + v); }, pwRs * pwRl / (pwRs + pwRl), -260, 260, 8192, -vb - 6, -vb + 30);
  }

  /* Output transformer primary + two 12" drivers, referred to one half primary. */
  void buildLoad() {
    const double k = (8000.0 / 4) / 16;            // 8k plate-to-plate onto a 16 ohm nominal load
    const double Rres = 5.5 * 16 * k, w0 = 2 * 3.141592653589793 * 88, Q = 3.5;
    const El els[] = {
      { 'C', 1, 0, 1.5e-9 },                              // winding capacitance
      { 'L', 1, 4, 40.0 / 4 }, { 'R', 4, 0, 50 },         // primary inductance + copper
      { 'R', 1, 0, 9e3 },                                 // plate resistance + core loss
      { 'R', 1, 2, 0.8 * 16 * k },                        // voice-coil resistance
      { 'L', 2, 3, 1.0e-3 * k }, { 'R', 2, 3, 64 * k },   // lossy coil inductance
      { 'R', 3, 0, Rres },                                // cone resonance
      { 'L', 3, 0, Rres / (Q * w0) },
      { 'C', 3, 0, Q / (Rres * w0) },
    };
    const In ins[] = { { 'I', 1 } };
    buildSS(load, 4, els, 10, ins, 1, 1, fo);
    double r = load.D[0]; for (int i = 0; i < load.ns; i++) r += load.C[i] * load.Q[i];
    Req = r;
    nSec = 2 * std::sqrt(16.0 / 8000);  // half-primary volts -> 16 ohm tap volts
  }

  static double cl(double x) { return std::min(0.9995, std::max(0.0005, x)); }
  static double par(double x, double y) { return x * y / (x + y); }

  void updateNetworks() {
    // loudness pots with 100p bright caps, 470k mixers (220p on channel I), Miller capacitance of V2a.
    // Channel II also feeds V2a through R1 100k + R17 470k; their junction goes off to the tremolo FET
    // through 0.22uF and the Intensity pot (220k) with R5 100k underneath.
    const double Rs1 = 38e3, Rs2 = 44e3, Rv = 1e6, Rm = 470e3;
    const double a1 = cl(audioTaper(p[VOL1])), a2 = cl(audioTaper(p[VOL2])), dI = cl(p[DEPTH]);
    const El mixEls[] = {
      { 'R', 4, 10, Rs1 }, { 'R', 10, 1, (1 - a1) * Rv }, { 'R', 1, 0, a1 * Rv }, { 'C', 10, 1, 100e-12 },
      { 'R', 5, 11, Rs2 }, { 'R', 11, 2, (1 - a2) * Rv }, { 'R', 2, 0, a2 * Rv }, { 'C', 11, 2, 100e-12 },
      { 'R', 1, 3, Rm }, { 'C', 1, 3, 220e-12 },      // R18, C6
      { 'R', 2, 3, Rm },                              // R16
      { 'R', 2, 6, 100e3 }, { 'R', 6, 3, Rm },        // R1, R17
      { 'C', 6, 7, 0.22e-6 },
      { 'R', 7, 8, (1 - dI) * 220e3 }, { 'R', 8, 9, dI * 220e3 }, { 'R', 9, 0, 100e3 },
      { 'C', 3, 0, 80e-12 },
    };
    const In mixIns[] = { { 'V', 4 }, { 'V', 5 }, { 'I', 8 } };
    buildSS(mix, 11, mixEls, 18, mixIns, 3, 3, fo);
    buildSS(mixD, 11, mixEls, 18, mixIns, 3, 8, fo);
    fetZm = mix.D[2]; fetZd = mixD.D[2];
    for (int i = 0; i < 5; i++) { fetZm += mix.C[i] * mix.Q[3 * i + 2]; fetZd += mixD.C[i] * mix.Q[3 * i + 2]; }
    const double z1 = par(a1 * Rv, Rs1 + (1 - a1) * Rv), z2 = par(a2 * Rv, Rs2 + (1 - a2) * Rv);
    const double r = par(par(Rm + z1, Rm + z2), Rm + 100e3 + z2);
    if (!haveClampV2 || std::fabs(r - rth) > 0.003 * rth) {
      rth = r; haveClampV2 = true;
      const LUT& t = tabIg2;
      buildClamp(clampV2, [&t](double v) { return lut(t, v); }, r, -150, 150, 8192, -4, 8);
    }
    // tone stack
    const double t = cl(p[TREBLE]), b = cl(audioTaper(p[BASS])), m = cl(p[MID]);
    const El toneEls[] = {
      { 'R', 7, 1, 1300 },
      { 'C', 1, 2, 250e-12 },
      { 'R', 2, 3, (1 - t) * 250e3 }, { 'R', 3, 4, t * 250e3 },
      { 'R', 1, 5, 56e3 },
      { 'C', 5, 4, 22e-9 }, { 'R', 4, 6, b * 1e6 },
      { 'C', 5, 6, 22e-9 }, { 'R', 6, 0, m * 25e3 },
    };
    const In vIn7[] = { { 'V', 7 } };
    buildSS(tone, 7, toneEls, 9, vIn7, 1, 3, fo);
    // negative feedback: 27k from the 16 ohm tap, 5k presence pot, 0.1uF on its wiper
    const double pr = cl(p[PRESENCE]);
    const El fbEls[] = {
      { 'R', 3, 1, 27e3 }, { 'R', 1, 2, (1 - pr) * 5e3 }, { 'R', 2, 0, pr * 5e3 },
      { 'C', 2, 0, 0.1e-6 },
    };
    const In vIn3[] = { { 'V', 3 } };
    buildSS(fb, 3, fbEls, 4, vIn3, 1, 1, fo);
    aF = (fb.D[0] + fb.C[0] * fb.Q[0]) * nSec;
    // tremolo oscillator V6a: plate -> 32n/(68k + Speed) -> 22n/503k -> 22n/1.47M -> grid,
    // 47n plate to cathode, 560R cathode (plus 2k2 when the footswitch stops it)
    const El oscEls[] = {
      { 'R', 1, 0, OSC::RL }, { 'C', 1, 5, OSC::Cpk }, { 'R', 5, 0, OSC::Rk + (tremOn ? 0 : OSC::Rfs) },
      { 'C', 1, 2, OSC::C1 }, { 'R', 2, 0, OSC::R1 + OSC::Rpot * audioTaper(1 - p[SPEED]) },
      { 'C', 2, 3, OSC::C2 }, { 'R', 3, 0, OSC::R2 },
      { 'C', 3, 4, OSC::C3 }, { 'R', 4, 0, OSC::R3 },
    };
    const In oscIns[] = { { 'I', 1 }, { 'I', 4 }, { 'I', 5 } };
    buildSS(osc, 5, oscEls, 9, oscIns, 3, 4, fs);
    buildSS(oscK, 5, oscEls, 9, oscIns, 3, 5, fs);
    buildSS(oscP, 5, oscEls, 9, oscIns, 3, 1, fs);
    double rz = osc.D[1]; for (int i = 0; i < 4; i++) rz += osc.C[i] * osc.Q[3 * i + 1];
    oscRz = rz;
    buildClamp(clampOsc, igT, rz, -300, 300, 8192, -3, 4);
    dirty = false;
  }

  /* Pick the transformer secondary voltage that puts ~445 V on the plates at idle. */
  void calibrateSupply() {
    const double iLoad = 2 * idle + 0.020;
    auto sim = [&](double vpk) {
      double vb = 445, sum = 0; int n = 0; const double dts = 1.0 / 48000;
      for (int i = 0; i < 24000; i++) {
        const double vac = vpk * std::fabs(std::sin(2 * 3.141592653589793 * 50 * i * dts));
        vb += dts / C1 * ((vac > vb ? (vac - vb) / Rrect : 0) - iLoad);
        if (i >= 14400) { sum += vb; n++; }
      }
      return sum / n;
    };
    double vpk = 480;
    for (int k = 0; k < 4; k++) vpk += 445 - sim(vpk);
    Vpk = vpk;
  }

  /* One base-rate step of the phase-shift oscillator; leaves the plate voltage in oscVp. */
  void stepLFO() {
    SS& o = osc; double* x = o.x; double* z = o.z; const double* P = o.P; const double* Q = o.Q; const double* C = o.C;
    const SS& ck = oscK; const SS& cp = oscP;
    const double ip = oscIp, u0 = OSC::Vs / OSC::RL - ip, a0 = o.up[0] + u0, a1 = o.up[1], a2 = o.up[2] + ip;
    double yg = o.D[0] * u0 + o.D[2] * ip, yk = ck.D[0] * u0 + ck.D[2] * ip;
    for (int i = 0; i < 4; i++) {
      double a = Q[3 * i] * a0 + Q[3 * i + 1] * a1 + Q[3 * i + 2] * a2; const int k = i * 4;
      for (int j = 0; j < 4; j++) a += P[k + j] * x[j];
      z[i] = a; yg += C[i] * a; yk += ck.C[i] * a;
    }
    const double vgk = lut(clampOsc, yg - yk), u1 = (vgk + yk - yg) / oscRz;   // grid conduction limits the swing
    double vp = cp.D[0] * u0 + cp.D[1] * u1 + cp.D[2] * ip;
    for (int i = 0; i < 4; i++) { x[i] = z[i] + Q[3 * i + 1] * u1; vp += cp.C[i] * x[i]; }
    o.up[0] = u0; o.up[1] = u1; o.up[2] = ip;
    oscIp = triIp(vgk, vp - yk, ECC83);
    oscVp = vp;
  }

  void smooth() {
    for (int k = 0; k < NUM_KNOBS; k++) {
      const double d = tgt[k] - p[k];
      if (d != 0) { p[k] = std::fabs(d) < 1e-4 ? tgt[k] : p[k] + d * 0.25; dirty = true; }
    }
    if (dirty) updateNetworks();
  }
};

}  // namespace bb
