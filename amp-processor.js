/* 1962-style 2x12 combo amplifier model — circuit-level DSP core + AudioWorklet wrapper.
 *
 * Signal path (all oversampled):
 *   input jacks -> V1 (ECC83, two halves, shared 820R/250uF cathode)
 *   -> two 1M loudness pots + 270k mixers (500pF on the High Treble side) + Miller C
 *   -> V2a gain stage DC-coupled to V2b cathode follower (solved together)
 *   -> FMV tone stack (250pF / 56k / 22n / 22n, 250k-1M-25k)
 *   -> long-tail phase inverter (82k/100k)
 *   Tremolo (from the 1962 reissue schematic): V6a phase-shift oscillator -> V6b cathode
 *   follower -> J174 FET shunting the R1/R17 junction of channel II's mixing network
 *   through 0.22uF and the Intensity pot. Channel II only.
 *   -> 0.1uF couplers with grid conduction (bias shift / blocking)
 *   -> push-pull KT66 into a reactive speaker load through the output transformer,
 *      27k/5k negative feedback + presence, solved implicitly every sample
 *   -> GZ34 rectifier + reservoir + choke-fed screens (sag and ripple)
 *
 * Linear sub-circuits are built from netlists (MNA -> state space -> trapezoidal).
 * Valve stages are solved offline on their load lines (Koren models) into tables.
 */
'use strict';

/* ------------------------------------------------------------------ tables */
function makeLUT(n, x0, x1, fn) {
  const d = new Float64Array(n + 1);
  const dx = (x1 - x0) / (n - 1);
  for (let i = 0; i < n; i++) d[i] = fn(x0 + i * dx);
  d[n] = d[n - 1];
  return { d, x0, x1, inv: 1 / dx, n };
}
function lut(t, x) {
  const p = (x - t.x0) * t.inv;
  if (p <= 0) return t.d[0];
  const m = t.n - 1;
  if (p >= m) return t.d[m];
  const i = p | 0, a = t.d[i];
  return a + (t.d[i + 1] - a) * (p - i);
}
let SL = 0; // slope returned alongside lutS()
function lutS(t, x) {
  const p = (x - t.x0) * t.inv;
  if (p <= 0) { SL = 0; return t.d[0]; }
  const m = t.n - 1;
  if (p >= m) { SL = 0; return t.d[m]; }
  const i = p | 0, a = t.d[i], s = t.d[i + 1] - a;
  SL = s * t.inv;
  return a + s * (p - i);
}

/* Inverse of m(v) = v + Rp*ig(v): the grid voltage reached when an open-circuit
 * voltage x drives a conducting grid through resistance Rp. Identity below vA. */
function buildClamp(ig, Rp, lo, hi, n, vA, vB, reuse) {
  const M = 2048;
  const mv = buildClamp._m || (buildClamp._m = new Float64Array(M + 1));
  const dv = (vB - vA) / M;
  for (let j = 0; j <= M; j++) { const v = vA + j * dv; mv[j] = v + Rp * ig(v); }
  const t = reuse || { d: new Float64Array(n + 1), x0: lo, x1: hi, inv: (n - 1) / (hi - lo), n };
  const dx = (hi - lo) / (n - 1), d = t.d;
  const m0 = mv[0], endSlope = dv / (mv[M] - mv[M - 1]);
  let j = 0;
  for (let i = 0; i < n; i++) {
    const x = lo + i * dx;
    if (x <= m0) { d[i] = x - (m0 - vA); continue; }
    while (j < M - 1 && mv[j + 1] < x) j++;
    if (x >= mv[M]) { d[i] = vB + (x - mv[M]) * endSlope; continue; }
    d[i] = vA + (j + (x - mv[j]) / (mv[j + 1] - mv[j])) * dv;
  }
  d[n] = d[n - 1];
  return t;
}

/* ------------------------------------------------------------ valve models */
const ECC83 = { mu: 100, ex: 1.4, kg1: 1060, kp: 600, kvb: 300 };
const KT66 = { mu: 8.2, ex: 1.35, kg1: 1350, kg2: 4500, kp: 45, kvb: 13 };

function triIp(vgk, vpk, T) {
  if (vpk <= 0) return 0;
  const t = T.kp * (1 / T.mu + vgk / Math.sqrt(T.kvb + vpk * vpk));
  const e1 = (vpk / T.kp) * (t > 30 ? t : Math.log1p(Math.exp(t)));
  return 2 * Math.pow(e1, T.ex) / T.kg1;
}
function gridI(vgk, g, s, vo) {
  const x = (vgk + vo) / s;
  const sp = x > 30 ? vgk + vo : s * Math.log1p(Math.exp(x));
  return g * sp * Math.sqrt(sp);
}
const igT = (v) => gridI(v, 5e-4, 0.12, 0.2);   // ECC83 grid
const igP = (v) => gridI(v, 1.5e-3, 0.3, 0.0);  // KT66 grid

function bisect(f, lo, hi, n) { // f decreasing; returns root
  for (let i = 0; i < n; i++) { const m = 0.5 * (lo + hi); if (f(m) > 0) lo = m; else hi = m; }
  return 0.5 * (lo + hi);
}
/* plate current of a triode with plate load Rp from supply Vs (cathode fixed) */
function loadLine(vgk, Vs, Rp, T) {
  if (Vs <= 0) return 0;
  return bisect((I) => triIp(vgk, Vs - I * Rp, T) - I, 0, Vs / Rp, 26);
}

/* --------------------------------------------- netlist -> discrete state space */
function solveLin(A, B, N, ne) { // Gauss-Jordan, A (N x N) and B (N x ne) as arrays of rows
  for (let c = 0; c < N; c++) {
    let p = c, best = Math.abs(A[c][c]);
    for (let r = c + 1; r < N; r++) { const v = Math.abs(A[r][c]); if (v > best) { best = v; p = r; } }
    if (p !== c) { let t = A[p]; A[p] = A[c]; A[c] = t; t = B[p]; B[p] = B[c]; B[c] = t; }
    const piv = 1 / A[c][c];
    for (let k = c; k < N; k++) A[c][k] *= piv;
    for (let k = 0; k < ne; k++) B[c][k] *= piv;
    for (let r = 0; r < N; r++) {
      if (r === c) continue;
      const f = A[r][c];
      if (f === 0) continue;
      for (let k = c; k < N; k++) A[r][k] -= f * A[c][k];
      for (let k = 0; k < ne; k++) B[r][k] -= f * B[c][k];
    }
  }
  return B;
}
function zeros(r, c) { const m = []; for (let i = 0; i < r; i++) m.push(new Float64Array(c)); return m; }

/* els: {t:'R'|'C'|'L', a, b, v}; ins: {t:'V'|'I', a}; node 0 is ground.
 * States: capacitor voltages (a-b) then inductor currents (a->b). */
function buildSS(nNodes, els, ins, outNode, fs, prev) {
  const caps = [], inds = [];
  for (const e of els) { if (e.t === 'C') caps.push(e); else if (e.t === 'L') inds.push(e); }
  const nc = caps.length, ns = nc + inds.length, ni = ins.length, ne = ns + ni;
  let nvs = nc; for (const e of ins) if (e.t === 'V') nvs++;
  const N = nNodes + nvs;
  const A = zeros(N, N), B = zeros(N, ne);
  for (const e of els) {
    if (e.t !== 'R') continue;
    const g = 1 / Math.max(e.v, 1), a = e.a - 1, b = e.b - 1;
    if (a >= 0) A[a][a] += g;
    if (b >= 0) A[b][b] += g;
    if (a >= 0 && b >= 0) { A[a][b] -= g; A[b][a] -= g; }
  }
  let row = nNodes;
  const vsrc = (a, b, exc) => {
    if (a >= 0) { A[row][a] = 1; A[a][row] = 1; }
    if (b >= 0) { A[row][b] = -1; A[b][row] = -1; }
    B[row][exc] = 1; row++;
  };
  caps.forEach((e, k) => vsrc(e.a - 1, e.b - 1, k));
  inds.forEach((e, k) => {
    if (e.a > 0) B[e.a - 1][nc + k] -= 1;
    if (e.b > 0) B[e.b - 1][nc + k] += 1;
  });
  ins.forEach((e, k) => {
    if (e.t === 'V') vsrc(e.a - 1, -1, ns + k); else B[e.a - 1][ns + k] += 1;
  });
  const X = solveLin(A, B, N, ne);
  const Mx = zeros(ns, ne);
  caps.forEach((e, k) => { for (let j = 0; j < ne; j++) Mx[k][j] = X[nNodes + k][j] / e.v; });
  inds.forEach((e, k) => {
    for (let j = 0; j < ne; j++) {
      Mx[nc + k][j] = ((e.a > 0 ? X[e.a - 1][j] : 0) - (e.b > 0 ? X[e.b - 1][j] : 0)) / e.v;
    }
  });
  const h2 = 0.5 / fs;
  const L = zeros(ns, ns), R = zeros(ns, ns + ni);
  for (let i = 0; i < ns; i++) {
    for (let j = 0; j < ns; j++) { L[i][j] = (i === j ? 1 : 0) - h2 * Mx[i][j]; R[i][j] = (i === j ? 1 : 0) + h2 * Mx[i][j]; }
    for (let j = 0; j < ni; j++) R[i][ns + j] = h2 * Mx[i][ns + j];
  }
  const S = solveLin(L, R, ns, ns + ni);
  const ss = prev || { ns, ni, P: new Float64Array(ns * ns), Q: new Float64Array(ns * ni), C: new Float64Array(ns), D: new Float64Array(ni), x: new Float64Array(ns), z: new Float64Array(ns), up: new Float64Array(ni) };
  for (let i = 0; i < ns; i++) {
    for (let j = 0; j < ns; j++) ss.P[i * ns + j] = S[i][j];
    for (let j = 0; j < ni; j++) ss.Q[i * ni + j] = S[i][ns + j];
    ss.C[i] = X[outNode - 1][i];
  }
  for (let j = 0; j < ni; j++) ss.D[j] = X[outNode - 1][ns + j];
  return ss;
}
/* single-input step: x <- P x + Q (u_prev + u), y = C x + D u */
function ssStep(s, u) {
  const ns = s.ns, P = s.P, Q = s.Q, x = s.x, z = s.z, uu = s.up[0] + u;
  let y = s.D[0] * u;
  for (let i = 0; i < ns; i++) {
    let a = Q[i] * uu; const o = i * ns;
    for (let j = 0; j < ns; j++) a += P[o + j] * x[j];
    z[i] = a; y += s.C[i] * a;
  }
  for (let i = 0; i < ns; i++) x[i] = z[i];
  s.up[0] = u;
  return y;
}
function ssStep2(s, u0, u1) {
  const ns = s.ns, P = s.P, Q = s.Q, x = s.x, z = s.z, a0 = s.up[0] + u0, a1 = s.up[1] + u1;
  let y = s.D[0] * u0 + s.D[1] * u1;
  for (let i = 0; i < ns; i++) {
    let a = Q[2 * i] * a0 + Q[2 * i + 1] * a1; const o = i * ns;
    for (let j = 0; j < ns; j++) a += P[o + j] * x[j];
    z[i] = a; y += s.C[i] * a;
  }
  for (let i = 0; i < ns; i++) x[i] = z[i];
  s.up[0] = u0; s.up[1] = u1;
  return y;
}

/* ------------------------------------------------------------- oversampling */
function kaiserSinc(taps, fc, beta) {
  const h = new Float64Array(taps), mid = (taps - 1) / 2;
  const i0 = (x) => { let s = 1, t = 1; for (let k = 1; k < 30; k++) { t *= (x / (2 * k)) * (x / (2 * k)); s += t; } return s; };
  let sum = 0;
  for (let i = 0; i < taps; i++) {
    const n = i - mid, r = n / mid;
    const w = i0(beta * Math.sqrt(Math.max(0, 1 - r * r))) / i0(beta);
    h[i] = (n === 0 ? 2 * fc : Math.sin(2 * Math.PI * fc * n) / (Math.PI * n)) * w;
    sum += h[i];
  }
  for (let i = 0; i < taps; i++) h[i] /= sum;
  return h;
}

/* Tremolo, values from the 1962 reissue schematic (V6 ECC83 + J174). */
const OSC = { Vs: 320, RL: 150e3, Rk: 560, Rfs: 2200, Cpk: 47e-9, C1: 32e-9, R1: 68e3, Rpot: 2e6, C2: 22e-9, R2: 503e3, C3: 22e-9, R3: 1.47e6 };
const FET = { vp: 7, idss: 0.045 };   // J174, P-channel, source grounded
function fetI(v, vg) {                 // current from the drain node into the FET
  const beta = FET.idss / (FET.vp * FET.vp);
  if (v <= 0) {
    const vov = FET.vp - vg, vds = -v;
    if (vov <= 0) return 0;
    return -(vds < vov ? beta * (2 * vov - vds) * vds : beta * vov * vov);
  }
  const vov = FET.vp - vg + v;
  if (vov <= 0) return 0;
  return v < vov ? beta * (2 * vov - v) * v : beta * vov * vov;
}

const audioTaper = (x) => (Math.pow(81, x) - 1) / 80; // 10% at half rotation

/* ================================================================ amp core */
class AmpCore {
  constructor(fs, os) {
    this.fs = fs; this.os = os; this.fo = fs * os; this.dt = 1 / this.fo;

    // oversampling filters
    const per = 20; this.per = per;
    this.h = kaiserSinc(per * os + 1, 0.47 / os, 7.9);
    this.hLen = this.h.length;
    this.upBuf = new Float64Array(64); this.upPos = 0;
    this.dnBuf = new Float64Array(1024); this.dnPos = 0;

    this.buildValves();

    // knob state (0..1) and targets
    this.p = { vol1: 0.6, vol2: 0.0, treble: 0.6, mid: 0.6, bass: 0.5, presence: 0.5, speed: 0.5, depth: 0 };
    this.tgt = Object.assign({}, this.p);
    this.mode = 0;        // 0 high treble, 1 normal, 2 jumpered
    this.sens = 1;        // 1 = high-sensitivity jack, 0.5 = low
    this.inScale = 2.0;   // volts at digital full scale
    this.dirty = true;

    // dynamic state
    this.vk1 = this.vk1Q; this.hp1 = 0; this.hp2 = 0;
    this.cPI = 0; this.cA = 0; this.cB = 0; this.vF = 0;
    this.Vh = 0; this.uL = 0; this.VhPrev = 0;
    this.Vb = 445; this.Vg2 = 438; this.acPh = 0;
    // tremolo valve state
    this.oscIp = 0.0005; this.oscVp = 200; this.tremOn = true; this.fetK = 0.035; this.fetV = 0; this.fetG = 20;
    this.sc = Math.pow(this.Vg2 / KT66.kp, KT66.ex) * 2 / KT66.kg1;
    this.hpC = 1 - Math.exp(-2 * Math.PI * 7 / this.fo);
    this.iters = 0;

    this.buildLoad();
    this.updateNetworks();
    this.calibrateSupply();
    // Let the tremolo oscillator build up at its slowest setting (smallest swing) and set the gate
    // trimmer (RV2) so the LFO peak still pinches the FET off there; then return to the knob setting.
    const spd = this.p.speed; this.p.speed = 0; this.updateNetworks();
    this.osc.x[0] = 60; for (let i = 0; i < fs * 5; i++) this.stepLFO();
    let top = 0; for (let i = 0; i < fs; i++) { this.stepLFO(); if (this.oscVp > top) top = this.oscVp; }
    this.fetK = (FET.vp + 0.6) / top;
    this.p.speed = spd; this.updateNetworks();
    for (let i = 0; i < fs * 2; i++) this.stepLFO();
    const z = new Float32Array(256), o = new Float32Array(256);
    for (let i = 0; i < Math.ceil(fs * 0.3 / 256); i++) this.process(z, o, 256);
  }

  buildValves() {
    const T = ECC83;
    /* ---- V1: both halves, 100k plates, shared 820R / 250uF cathode ---- */
    const Vb1 = 265, Rp = 100e3, Rk = 820;
    this.Rp1 = Rp; this.Rk1 = Rk; this.Ck1 = 250e-6;
    this.vk1Q = bisect((vk) => 2 * loadLine(-vk, Vb1 - vk, Rp, T) - vk / Rk, 0, 6, 30);
    this.tabV1 = makeLUT(4096, -9, 4, (v) => loadLine(v, Vb1 - this.vk1Q, Rp, T));
    this.ip1Q = lut(this.tabV1, -this.vk1Q);
    this.clampIn = buildClamp(igT, 68e3, -20, 20, 8192, -3, 4);

    /* ---- V2a (100k / 820R unbypassed) DC-coupled to V2b follower (100k) ---- */
    const Vb2 = 300, Rk2 = 820, Rcf = 100e3;
    const solve = (vg) => {
      let i1 = 0, vk2 = 0;
      const vp = bisect((Vp) => {
        i1 = bisect((I) => triIp(vg - Rk2 * I, Vp - Rk2 * I, T) - I, 0, 0.012, 24);
        vk2 = bisect((Vk) => triIp(Vp - Vk, Vb2 - Vk, T) + igT(Vp - Vk) - Vk / Rcf, 0, Vb2, 24);
        return (Vb2 - Vp) / Rp - i1 - igT(Vp - vk2);
      }, 0, Vb2, 24);
      void vp;
      return [vk2, igT(vg - Rk2 * i1)];
    };
    const n2 = 2048, cf = new Float64Array(n2), ig2 = new Float64Array(n2);
    for (let i = 0; i < n2; i++) { const r = solve(-10 + i * 18 / (n2 - 1)); cf[i] = r[0]; ig2[i] = r[1]; }
    const q2 = solve(0)[0];
    this.tabV2 = makeLUT(n2, -10, 8, () => 0); this.tabIg2 = makeLUT(n2, -10, 8, () => 0);
    for (let i = 0; i < n2; i++) { this.tabV2.d[i] = cf[i] - q2; this.tabIg2.d[i] = ig2[i]; }
    this.tabV2.d[n2] = this.tabV2.d[n2 - 1]; this.tabIg2.d[n2] = this.tabIg2.d[n2 - 1];
    this.clampV2 = null;

    /* ---- long-tail phase inverter: 82k / 100k, 470R bias, 10k + 5k tail ---- */
    const Vb3 = 320, Ra = 82e3, Rb = 100e3, Rbias = 470, Rtail = 15e3, Rkk = Rbias + Rtail;
    const itQ = bisect((It) => loadLine(-Rbias * It, Vb3 - Rkk * It, Ra, T) + loadLine(-Rbias * It, Vb3 - Rkk * It, Rb, T) - It, 0, 0.006, 30);
    const vxQ = itQ * Rtail;
    const piSolve = (vd) => {
      let ia = 0, ib = 0, ig = 0;
      bisect((Vk) => {
        ia = loadLine(vxQ + vd - Vk, Vb3 - Vk, Ra, T);
        ib = loadLine(vxQ - Vk, Vb3 - Vk, Rb, T);
        ig = igT(vxQ + vd - Vk);
        return ia + ib + ig - Vk / Rkk;
      }, 0, Vb3, 26);
      return [ia, ib, ig];
    };
    const n3 = 1024, q3 = piSolve(0);
    this.tabPIa = makeLUT(n3, -30, 30, () => 0); this.tabPIb = makeLUT(n3, -30, 30, () => 0);
    const tig = makeLUT(n3, -30, 30, () => 0);
    for (let i = 0; i < n3; i++) {
      const r = piSolve(-30 + i * 60 / (n3 - 1));
      this.tabPIa.d[i] = -Ra * (r[0] - q3[0]); this.tabPIb.d[i] = -Rb * (r[1] - q3[1]); tig.d[i] = r[2];
    }
    this.tabPIa.d[n3] = this.tabPIa.d[n3 - 1]; this.tabPIb.d[n3] = this.tabPIb.d[n3 - 1]; tig.d[n3] = tig.d[n3 - 1];
    // tone stack -> 22n -> PI grid, 1M leak
    this.piRs = 120e3; this.piRl = 1e6; this.piC = 22e-9;
    this.piK = this.piRl / (this.piRs + this.piRl);
    this.clampPI = buildClamp((v) => lut(tig, v), this.piRs * this.piRl / (this.piRs + this.piRl), -150, 150, 8192, -30, 30);

    /* ---- KT66 pair ---- */
    const P = KT66;
    this.tabL = makeLUT(4096, -0.5, 0.08, (x) => {
      const t = P.kp * (1 / P.mu + x);
      return Math.pow(t > 30 ? t : Math.log1p(Math.exp(t)), P.ex);
    });
    const idle = 0.036, VbN = 445, Vg2N = 438;
    const scN = Math.pow(Vg2N / P.kp, P.ex) * 2 / P.kg1;
    this.vBias = bisect((vb) => idle - scN * lut(this.tabL, vb / Vg2N) * Math.atan(VbN / P.kvb), -90, -10, 40);
    this.idle = idle;
    this.pwRs = 50e3; this.pwRl = 220e3; this.pwC = 0.1e-6;
    this.pwK = this.pwRl / (this.pwRs + this.pwRl);
    const vb = this.vBias;
    this.clampPw = buildClamp((v) => igP(vb + v), this.pwRs * this.pwRl / (this.pwRs + this.pwRl), -260, 260, 8192, -vb - 6, -vb + 30);
  }

  /* Output transformer primary + two 12" drivers, referred to one half primary. */
  buildLoad() {
    const k = (8000 / 4) / 16;            // 8k plate-to-plate onto a 16 ohm nominal load
    const Rres = 5.5 * 16 * k, w0 = 2 * Math.PI * 88, Q = 3.5;
    this.load = buildSS(4, [
      { t: 'C', a: 1, b: 0, v: 1.5e-9 },                       // winding capacitance
      { t: 'L', a: 1, b: 4, v: 40 / 4 }, { t: 'R', a: 4, b: 0, v: 50 }, // primary inductance + copper
      { t: 'R', a: 1, b: 0, v: 9e3 },                          // plate resistance + core loss
      { t: 'R', a: 1, b: 2, v: 0.8 * 16 * k },                 // voice-coil resistance
      { t: 'L', a: 2, b: 3, v: 1.0e-3 * k }, { t: 'R', a: 2, b: 3, v: 64 * k }, // lossy coil inductance
      { t: 'R', a: 3, b: 0, v: Rres },                         // cone resonance
      { t: 'L', a: 3, b: 0, v: Rres / (Q * w0) },
      { t: 'C', a: 3, b: 0, v: Q / (Rres * w0) },
    ], [{ t: 'I', a: 1 }], 1, this.fo);
    const s = this.load;
    let r = s.D[0]; for (let i = 0; i < s.ns; i++) r += s.C[i] * s.Q[i];
    this.Req = r;
    this.nSec = 2 * Math.sqrt(16 / 8000);  // half-primary volts -> 16 ohm tap volts
  }

  updateNetworks() {
    const p = this.p, fo = this.fo, cl = (x) => Math.min(0.9995, Math.max(0.0005, x));
    // loudness pots with 100p bright caps, 470k mixers (220p on channel I), Miller capacitance of V2a.
    // Channel II also feeds V2a through R1 100k + R17 470k; their junction goes off to the tremolo FET
    // through 0.22uF and the Intensity pot (220k) with R5 100k underneath.
    const Rs1 = 38e3, Rs2 = 44e3, Rv = 1e6, Rm = 470e3;
    const a1 = cl(audioTaper(p.vol1)), a2 = cl(audioTaper(p.vol2)), dI = cl(p.depth);
    const mixEls = [
      { t: 'R', a: 4, b: 10, v: Rs1 }, { t: 'R', a: 10, b: 1, v: (1 - a1) * Rv }, { t: 'R', a: 1, b: 0, v: a1 * Rv }, { t: 'C', a: 10, b: 1, v: 100e-12 },
      { t: 'R', a: 5, b: 11, v: Rs2 }, { t: 'R', a: 11, b: 2, v: (1 - a2) * Rv }, { t: 'R', a: 2, b: 0, v: a2 * Rv }, { t: 'C', a: 11, b: 2, v: 100e-12 },
      { t: 'R', a: 1, b: 3, v: Rm }, { t: 'C', a: 1, b: 3, v: 220e-12 },      // R18, C6
      { t: 'R', a: 2, b: 3, v: Rm },                                            // R16
      { t: 'R', a: 2, b: 6, v: 100e3 }, { t: 'R', a: 6, b: 3, v: Rm },         // R1, R17
      { t: 'C', a: 6, b: 7, v: 0.22e-6 },
      { t: 'R', a: 7, b: 8, v: (1 - dI) * 220e3 }, { t: 'R', a: 8, b: 9, v: dI * 220e3 }, { t: 'R', a: 9, b: 0, v: 100e3 },
      { t: 'C', a: 3, b: 0, v: 80e-12 },
    ];
    const mixIns = [{ t: 'V', a: 4 }, { t: 'V', a: 5 }, { t: 'I', a: 8 }];
    this.mix = buildSS(11, mixEls, mixIns, 3, fo, this.mix);
    this.mixD = buildSS(11, mixEls, mixIns, 8, fo, this.mixD);
    this.fetZm = this.mix.D[2]; this.fetZd = this.mixD.D[2];
    for (let i = 0; i < 5; i++) { this.fetZm += this.mix.C[i] * this.mix.Q[3 * i + 2]; this.fetZd += this.mixD.C[i] * this.mix.Q[3 * i + 2]; }
    const par = (x, y) => x * y / (x + y);
    const z1 = par(a1 * Rv, Rs1 + (1 - a1) * Rv), z2 = par(a2 * Rv, Rs2 + (1 - a2) * Rv);
    const rth = par(par(Rm + z1, Rm + z2), Rm + 100e3 + z2);
    if (!this.clampV2 || Math.abs(rth - this.rth) > 0.003 * this.rth) {
      this.rth = rth;
      const t = this.tabIg2;
      this.clampV2 = buildClamp((v) => lut(t, v), rth, -150, 150, 8192, -4, 8, this.clampV2);
    }
    // tone stack
    const t = cl(p.treble), b = cl(audioTaper(p.bass)), m = cl(p.mid);
    this.tone = buildSS(7, [
      { t: 'R', a: 7, b: 1, v: 1300 },
      { t: 'C', a: 1, b: 2, v: 250e-12 },
      { t: 'R', a: 2, b: 3, v: (1 - t) * 250e3 }, { t: 'R', a: 3, b: 4, v: t * 250e3 },
      { t: 'R', a: 1, b: 5, v: 56e3 },
      { t: 'C', a: 5, b: 4, v: 22e-9 }, { t: 'R', a: 4, b: 6, v: b * 1e6 },
      { t: 'C', a: 5, b: 6, v: 22e-9 }, { t: 'R', a: 6, b: 0, v: m * 25e3 },
    ], [{ t: 'V', a: 7 }], 3, fo, this.tone);
    // negative feedback: 27k from the 16 ohm tap, 5k presence pot, 0.1uF on its wiper
    const pr = cl(p.presence);
    this.fb = buildSS(3, [
      { t: 'R', a: 3, b: 1, v: 27e3 }, { t: 'R', a: 1, b: 2, v: (1 - pr) * 5e3 }, { t: 'R', a: 2, b: 0, v: pr * 5e3 },
      { t: 'C', a: 2, b: 0, v: 0.1e-6 },
    ], [{ t: 'V', a: 3 }], 1, fo, this.fb);
    this.aF = (this.fb.D[0] + this.fb.C[0] * this.fb.Q[0]) * this.nSec;
    // tremolo oscillator V6a: plate -> 32n/(68k + Speed) -> 22n/503k -> 22n/1.47M -> grid,
    // 47n plate to cathode, 560R cathode (plus 2k2 when the footswitch stops it)
    const oscEls = [
      { t: 'R', a: 1, b: 0, v: OSC.RL }, { t: 'C', a: 1, b: 5, v: OSC.Cpk }, { t: 'R', a: 5, b: 0, v: OSC.Rk + (this.tremOn ? 0 : OSC.Rfs) },
      { t: 'C', a: 1, b: 2, v: OSC.C1 }, { t: 'R', a: 2, b: 0, v: OSC.R1 + OSC.Rpot * audioTaper(1 - p.speed) },
      { t: 'C', a: 2, b: 3, v: OSC.C2 }, { t: 'R', a: 3, b: 0, v: OSC.R2 },
      { t: 'C', a: 3, b: 4, v: OSC.C3 }, { t: 'R', a: 4, b: 0, v: OSC.R3 },
    ], oscIns = [{ t: 'I', a: 1 }, { t: 'I', a: 4 }, { t: 'I', a: 5 }];
    this.osc = buildSS(5, oscEls, oscIns, 4, this.fs, this.osc);
    this.oscK = buildSS(5, oscEls, oscIns, 5, this.fs, this.oscK);
    this.oscP = buildSS(5, oscEls, oscIns, 1, this.fs, this.oscP);
    let rz = this.osc.D[1]; for (let i = 0; i < 4; i++) rz += this.osc.C[i] * this.osc.Q[3 * i + 1];
    this.oscRz = rz;
    this.clampOsc = buildClamp(igT, rz, -300, 300, 8192, -3, 4, this.clampOsc);
    this.dirty = false;
  }

  /* Pick the transformer secondary voltage that puts ~445 V on the plates at idle. */
  calibrateSupply() {
    this.Rrect = 90; this.C1 = 32e-6; this.Rch = 250; this.C2 = 32e-6; this.iPre = 0.006;
    const iLoad = 2 * this.idle + 0.020;
    const sim = (vpk) => {
      let vb = 445, sum = 0, n = 0; const dt = 1 / 48000;
      for (let i = 0; i < 24000; i++) {
        const vac = vpk * Math.abs(Math.sin(2 * Math.PI * 50 * i * dt));
        vb += dt / this.C1 * ((vac > vb ? (vac - vb) / this.Rrect : 0) - iLoad);
        if (i >= 14400) { sum += vb; n++; }
      }
      return sum / n;
    };
    let vpk = 480;
    for (let k = 0; k < 4; k++) vpk += 445 - sim(vpk);
    this.Vpk = vpk;
  }

  /* One base-rate step of the phase-shift oscillator; leaves the plate voltage in this.oscVp. */
  stepLFO() {
    const o = this.osc, x = o.x, z = o.z, P = o.P, Q = o.Q, C = o.C, ck = this.oscK, cp = this.oscP;
    const ip = this.oscIp, u0 = OSC.Vs / OSC.RL - ip, a0 = o.up[0] + u0, a1 = o.up[1], a2 = o.up[2] + ip;
    let yg = o.D[0] * u0 + o.D[2] * ip, yk = ck.D[0] * u0 + ck.D[2] * ip;
    for (let i = 0; i < 4; i++) {
      let a = Q[3 * i] * a0 + Q[3 * i + 1] * a1 + Q[3 * i + 2] * a2; const k = i * 4;
      for (let j = 0; j < 4; j++) a += P[k + j] * x[j];
      z[i] = a; yg += C[i] * a; yk += ck.C[i] * a;
    }
    const vgk = lut(this.clampOsc, yg - yk), u1 = (vgk + yk - yg) / this.oscRz;   // grid conduction limits the swing
    let vp = cp.D[0] * u0 + cp.D[1] * u1 + cp.D[2] * ip;
    for (let i = 0; i < 4; i++) { x[i] = z[i] + Q[3 * i + 1] * u1; vp += cp.C[i] * x[i]; }
    o.up[0] = u0; o.up[1] = u1; o.up[2] = ip;
    this.oscIp = triIp(vgk, vp - yk, ECC83);
    this.oscVp = vp;
  }

  setKnob(name, v) { if (name in this.tgt) this.tgt[name] = Math.min(1, Math.max(0, v)); }

  smooth() {
    const p = this.p, t = this.tgt;
    for (const k in p) {
      const d = t[k] - p[k];
      if (d !== 0) { p[k] = Math.abs(d) < 1e-4 ? t[k] : p[k] + d * 0.25; this.dirty = true; }
    }
    if (this.dirty) this.updateNetworks();
  }

  /* inp/out: Float32Array at the base rate. Output is normalised half-primary voltage. */
  process(inp, out, n) {
    this.smooth();
    const os = this.os, per = this.per, h = this.h, hLen = this.hLen;
    const upBuf = this.upBuf, dnBuf = this.dnBuf;
    const dt = this.dt, P = KT66, ikvb = 1 / P.kvb, kgr = P.kg1 / P.kg2;
    const tabV1 = this.tabV1, clampIn = this.clampIn, clampV2 = this.clampV2, tabV2 = this.tabV2;
    const clampPI = this.clampPI, tabPIa = this.tabPIa, tabPIb = this.tabPIb, clampPw = this.clampPw, tabL = this.tabL;
    const mix = this.mix, tone = this.tone, fb = this.fb, load = this.load;
    const Rp1 = this.Rp1, ipQ = this.ip1Q, kCk = dt / this.Ck1, iRk = 1 / this.Rk1, hpC = this.hpC;
    const gB = this.mode === 1 ? 0 : this.sens * this.inScale, gN = this.mode === 0 ? 0 : this.sens * this.inScale;
    const piK = this.piK, piRs = 1 / this.piRs, piC = dt / this.piC;
    const pwK = this.pwK, pwRs = 1 / this.pwRs, pwC = dt / this.pwC;
    const aF = this.aF, Req = this.Req, vBias = this.vBias, nSec = this.nSec;
    const mx = mix.x, mz = mix.z, mP = mix.P, mQ = mix.Q, mC = mix.C, mD = mix.D, mup = mix.up, dC = this.mixD.C, dD = this.mixD.D;
    const Zm = this.fetZm, Zd = this.fetZd, fetK = this.fetK;
    let fetV = this.fetV, fetG = this.fetG;
    const acInc = 2 * Math.PI * 50 * dt, Vpk = this.Vpk, iRr = 1 / this.Rrect, iRch = 1 / this.Rch, kC1 = dt / this.C1, kC2 = dt / this.C2, iPre = this.iPre;
    const lx = load.x, lz = load.z, lP = load.P, lQ = load.Q, lC = load.C, lns = load.ns;
    const fx = fb.x, fP = fb.P[0], fQ = fb.Q[0], fC = fb.C[0];
    let vk1 = this.vk1, hp1 = this.hp1, hp2 = this.hp2, cPI = this.cPI, cA = this.cA, cB = this.cB, vF = this.vF;
    let Vh = this.Vh, uL = this.uL, VhPrev = this.VhPrev, Vb = this.Vb, Vg2 = this.Vg2, sc = this.sc;
    let acPh = this.acPh, upPos = this.upPos, dnPos = this.dnPos, iters = 0;

    for (let s = 0; s < n; s++) {
      upPos = (upPos + 1) & 63; upBuf[upPos] = inp[s];
      this.stepLFO();
      fetG = fetK * (this.oscVp - 1);   // V6b follower -> R40 / RV2 divider -> J174 gate
      for (let ph = 0; ph < os; ph++) {
        // ---- polyphase interpolation
        let x = 0;
        for (let k = 0, j = ph; k < per; k++, j += os) x += h[j] * upBuf[(upPos - k) & 63];
        x *= os;

        // ---- V1
        const ia = lut(tabV1, lut(clampIn, x * gB - vk1));
        const ib = lut(tabV1, lut(clampIn, x * gN - vk1));
        vk1 += kCk * (ia + ib - vk1 * iRk);
        let u1 = -Rp1 * (ia - ipQ), u2 = -Rp1 * (ib - ipQ);
        hp1 += hpC * (u1 - hp1); u1 -= hp1;
        hp2 += hpC * (u2 - hp2); u2 -= hp2;

        // ---- loudness pots + mixers, V2a + follower, tone stack
        // (the tremolo FET is solved implicitly against the node it shunts)
        const b0 = mup[0] + u1, b1 = mup[1] + u2, b2 = mup[2];
        let ym = mD[0] * u1 + mD[1] * u2, yd = dD[0] * u1 + dD[1] * u2;
        for (let i = 0; i < 5; i++) {
          const k = 5 * i, q = 3 * i;
          const a = mQ[q] * b0 + mQ[q + 1] * b1 + mQ[q + 2] * b2 + mP[k] * mx[0] + mP[k + 1] * mx[1] + mP[k + 2] * mx[2] + mP[k + 3] * mx[3] + mP[k + 4] * mx[4];
          mz[i] = a; ym += mC[i] * a; yd += dC[i] * a;
        }
        for (let it = 0; it < 6; it++) {
          const i0 = fetI(fetV, fetG), di = (fetI(fetV + 1e-3, fetG) - i0) * 1e3;
          const stp = (fetV - yd + Zd * i0) / (1 + Zd * di);
          fetV -= stp;
          if (stp < 1e-6 && stp > -1e-6) break;
        }
        const u3 = (fetV - yd) / Zd;
        for (let i = 0; i < 5; i++) mx[i] = mz[i] + mQ[3 * i + 2] * u3;
        mup[0] = u1; mup[1] = u2; mup[2] = u3;
        const vm = ym + Zm * u3;
        const v2 = lut(tabV2, lut(clampV2, vm));
        const vt = ssStep(tone, v2);

        // ---- phase inverter input coupler (grid conduction charges the 22n)
        const vn = lut(clampPI, piK * (vt - cPI) - vF) + vF;
        cPI += piC * (vt - cPI - vn) * piRs;

        // ---- power stage: implicit solve for the half-primary voltage
        let Veq = 0;
        for (let i = 0; i < lns; i++) {
          let a = lQ[i] * uL; const o = i * lns;
          for (let j = 0; j < lns; j++) a += lP[o + j] * lx[j];
          lz[i] = a; Veq += lC[i] * a;
        }
        const zF = fP * fx[0] + fQ * nSec * VhPrev, hF = fC * zF;
        let lo = -1e5, hi = 1e5, va = 0, vb = 0, gA = 0, gB2 = 0, IA = 0, IB = 0, LA = 0, LB = 0;
        const iVg2 = 1 / Vg2;
        for (let it = 0; it < 14; it++) {
          iters++;
          vF = aF * Vh + hF;
          const vd = vn - vF;
          va = lutS(tabPIa, vd); const sa = SL;
          vb = lutS(tabPIb, vd); const sb = SL;
          gA = lutS(clampPw, pwK * (vb - cA)); const sgA = SL * pwK * sb;
          gB2 = lutS(clampPw, pwK * (va - cB)); const sgB = SL * pwK * sa;
          LA = lutS(tabL, (vBias + gA) * iVg2); const sLA = SL * iVg2 * sgA;
          LB = lutS(tabL, (vBias + gB2) * iVg2); const sLB = SL * iVg2 * sgB;
          const pA = (Vb - Vh) * ikvb, pB = (Vb + Vh) * ikvb;
          const atA = pA > 0 ? Math.atan(pA) : 0, atB = pB > 0 ? Math.atan(pB) : 0;
          const dA = pA > 0 ? ikvb / (1 + pA * pA) : 0, dB = pB > 0 ? ikvb / (1 + pB * pB) : 0;
          IA = sc * LA * atA; IB = sc * LB * atB;
          const g = Vh - Req * (IA - IB) - Veq;
          const dg = 1 - Req * sc * ((-sLA * aF * atA - LA * dA) - (-sLB * aF * atB + LB * dB));
          if (g > 0) hi = Vh; else lo = Vh;
          let nx = Vh - g / dg;
          if (!(nx > lo && nx < hi)) nx = 0.5 * (lo + hi);
          const e = nx - Vh; Vh = nx;
          if (e < 2e-3 && e > -2e-3) break;
        }
        vF = aF * Vh + hF;
        const u = (Vh - Veq) / Req;
        for (let i = 0; i < lns; i++) lx[i] = lz[i] + lQ[i] * u;
        uL = u;
        fx[0] = zF + fQ * nSec * Vh; VhPrev = Vh;
        cA += pwC * (vb - cA - gA) * pwRs;
        cB += pwC * (va - cB - gB2) * pwRs;

        // ---- GZ34 + reservoir, choke-fed screens
        acPh += acInc; if (acPh > 6.283185307) acPh -= 6.283185307;
        let vac = Vpk * Math.sin(acPh); if (vac < 0) vac = -vac;
        const ich = (Vb - Vg2) * iRch;
        Vb += kC1 * ((vac > Vb ? (vac - Vb) * iRr : 0) - IA - IB - ich);
        Vg2 += kC2 * (ich - sc * (LA + LB) * kgr - iPre);
        sc = Math.pow(Vg2 / P.kp, P.ex) * 2 / P.kg1;

        dnPos = (dnPos + 1) & 1023; dnBuf[dnPos] = Vh;
      }
      // ---- decimation
      let y = 0;
      for (let j = 0; j < hLen; j++) y += h[j] * dnBuf[(dnPos - j) & 1023];
      out[s] = y * (1 / 430);
    }
    if (!(Vh === Vh) || !(vk1 === vk1) || !(fetV === fetV)) { this.resetState(); return; }
    this.vk1 = vk1; this.hp1 = hp1; this.hp2 = hp2; this.cPI = cPI; this.cA = cA; this.cB = cB; this.vF = vF;
    this.Vh = Vh; this.uL = uL; this.VhPrev = VhPrev; this.Vb = Vb; this.Vg2 = Vg2; this.sc = sc;
    this.acPh = acPh; this.fetV = fetV; this.fetG = fetG; this.upPos = upPos; this.dnPos = dnPos; this.iters = iters;
  }

  resetState() {
    this.vk1 = this.vk1Q; this.hp1 = this.hp2 = 0; this.cPI = this.cA = this.cB = 0; this.vF = 0;
    this.Vh = 0; this.uL = 0; this.VhPrev = 0; this.Vb = 445; this.Vg2 = 438;
    this.fetV = 0;
    this.sc = Math.pow(this.Vg2 / KT66.kp, KT66.ex) * 2 / KT66.kg1;
    for (const s of [this.mix, this.tone, this.fb, this.load]) { s.x.fill(0); s.up.fill(0); }
    this.upBuf.fill(0); this.dnBuf.fill(0);
  }
}

/* ============================================================ worklet glue */
if (typeof registerProcessor === 'function') {
  class AmpProcessor extends AudioWorkletProcessor {
    constructor(options) {
      super();
      const o = (options && options.processorOptions) || {};
      this.core = new AmpCore(sampleRate, o.os || 4);
      this.standby = false;
      this.scratch = new Float32Array(128);
      this.count = 0; this.peakIn = 0; this.peakOut = 0;
      this.port.onmessage = (e) => {
        const m = e.data, c = this.core;
        if (m.knobs) for (const k in m.knobs) c.setKnob(k, m.knobs[k]);
        if (m.mode !== undefined) c.mode = m.mode;
        if (m.sens !== undefined) c.sens = m.sens;
        if (m.inScale !== undefined) c.inScale = m.inScale;
        if (m.standby !== undefined) this.standby = m.standby;
        if (m.trem !== undefined && m.trem !== c.tremOn) { c.tremOn = m.trem; c.dirty = true; }
      };
      this.port.postMessage({ ready: true });
    }
    process(inputs, outputs) {
      const out = outputs[0][0], inp = inputs[0] && inputs[0][0];
      if (!out) return true;
      const n = out.length;
      if (this.scratch.length < n) this.scratch = new Float32Array(n);
      const x = this.scratch;
      if (inp && !this.standby) x.set(inp); else x.fill(0, 0, n);
      this.core.process(x, out, n);
      let pi = 0, po = 0;
      for (let i = 0; i < n; i++) { const a = Math.abs(x[i]), b = Math.abs(out[i]); if (a > pi) pi = a; if (b > po) po = b; }
      if (pi > this.peakIn) this.peakIn = pi;
      if (po > this.peakOut) this.peakOut = po;
      if (++this.count >= 16) {
        this.port.postMessage({ pin: this.peakIn, pout: this.peakOut, vb: this.core.Vb, it: this.core.iters / (n * this.core.os) });
        this.count = 0; this.peakIn = 0; this.peakOut = 0;
      }
      return true;
    }
  }
  registerProcessor('amp-1962', AmpProcessor);
}
if (typeof window !== 'undefined') window.AmpCore = AmpCore;
