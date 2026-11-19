// qPCR analysis shared by the run viewer (inlined into each run page) and the site index
// (loaded as analysis.js). Pure functions only: no DOM. Keep this the single copy of the
// Cq / Tm logic so the index export and the viewer always agree.
//
// run = the JSON written by qpcr_viewer.py: {plate, cycles, melt_temps, channels, wells, ...}
"use strict";
const QPCR = (() => {
  const ROWS = "ABCDEFGHIJKLMNOP";
  const wellsOf = run => {
    const out = [];
    for (let r = 0; r < run.plate.rows; r++) for (let c = 0; c < run.plate.cols; c++) out.push(ROWS[r] + (c + 1));
    return out;
  };
  const mean = a => a.length ? a.reduce((s, v) => s + v, 0) / a.length : NaN;
  const sd = a => { if (a.length < 2) return NaN; const m = mean(a); return Math.sqrt(a.reduce((s, v) => s + (v - m) ** 2, 0) / (a.length - 1)); };
  const median = a => { if (!a.length) return NaN; const s = [...a].sort((x, y) => x - y); const m = s.length >> 1; return s.length % 2 ? s[m] : (s[m - 1] + s[m]) / 2; };

  const nonzero = arr => !!arr && arr.some(v => v != null && v !== 0);
  const hasData = (ch, w) => nonzero(ch.amp[w]) || nonzero(ch.melt[w]);

  // Default baseline window: cycles 3-15, clipped for short runs.
  const defaultBaseline = cyc => [Math.min(3, cyc.length || 3), Math.min(15, Math.max(4, cyc.length - 5))];

  function linfit(xs, ys) {
    const n = xs.length;
    if (n < 2) return { a: n ? ys[0] : 0, b: 0 };
    let sx = 0, sy = 0, sxx = 0, sxy = 0;
    for (let i = 0; i < n; i++) { sx += xs[i]; sy += ys[i]; sxx += xs[i] * xs[i]; sxy += xs[i] * ys[i]; }
    const den = n * sxx - sx * sx;
    const b = den ? (n * sxy - sx * sy) / den : 0;
    return { a: (sy - b * sx) / n, b };
  }
  function findCq(corr, thr, cyc) {
    // last point below threshold; Cq interpolated between it and the next point
    let j = -1;
    for (let i = 0; i < corr.length; i++) if (corr[i] != null && corr[i] < thr) j = i;
    if (j < 0 || j >= corr.length - 1) return null;
    let k = j + 1; while (k < corr.length && corr[k] == null) k++;
    if (k >= corr.length) return null;
    const y0 = corr[j], y1 = corr[k];
    return cyc[j] + (thr - y0) / (y1 - y0) * (cyc[k] - cyc[j]);
  }
  function negDeriv(ys, temps) {
    // 5-point Savitzky-Golay first derivative (window shrinks at the ends)
    const n = ys.length, out = new Array(n).fill(null);
    if (n < 3) return out;
    for (let i = 0; i < n; i++) {
      const m = Math.min(2, i, n - 1 - i);
      if (m === 0) continue;
      let num = 0, den = 0, ok = true;
      for (let k = -m; k <= m; k++) {
        const v = ys[i + k]; if (v == null) { ok = false; break; }
        num += k * v; den += k * (temps[i + k] - temps[i]);
      }
      if (ok && den) out[i] = -num / den;
    }
    out[0] = out[1]; out[n - 1] = out[n - 2];
    return out;
  }
  function findPeaks(d, thr, temps) {
    const peaks = [];
    for (let i = 2; i < d.length - 2; i++) {
      const a = d[i - 1], b = d[i], c = d[i + 1];
      if (a == null || b == null || c == null) continue;
      if (b >= thr && b >= a && b > c) {
        // prominence: how far the peak rises above the higher of the two valleys
        // separating it from taller signal (or the curve ends)
        let lmin = b, rmin = b;
        for (let j = i - 1; j >= 0 && d[j] != null && d[j] <= b; j--) lmin = Math.min(lmin, d[j]);
        for (let j = i + 1; j < d.length && d[j] != null && d[j] <= b; j++) rmin = Math.min(rmin, d[j]);
        if (b - Math.max(lmin, rmin) < Math.max(0.5 * thr, 0.1 * b)) continue;
        const den = a - 2 * b + c;
        const off = den ? 0.5 * (a - c) / den : 0;
        const step = temps[i + 1] - temps[i];
        peaks.push({ t: temps[i] + off * step, h: b - 0.25 * (a - c) * off });
      }
    }
    peaks.sort((p, q) => q.h - p.h);
    const kept = [];
    for (const p of peaks) if (kept.every(k => Math.abs(k.t - p.t) >= 2)) kept.push(p);
    return kept;
  }

  // Analyse one channel of a run.
  // opts: {base0, base1} baseline cycles (default 3-15); thr / peakThr: null or omitted = auto.
  // Returns {res: {well: {cq, tm, peaks, corr, deriv, fit}}, thr, autoThr, pthr, autoPeak, base0, base1}.
  function analyze(run, chanIdx = 0, opts = {}) {
    const ch = run.channels[chanIdx] || { amp: {}, melt: {} };
    const CYC = run.cycles || [], TEMPS = run.melt_temps || [];
    const [d0, d1] = defaultBaseline(CYC);
    const base0 = opts.base0 ?? d0, base1 = opts.base1 ?? d1;
    const res = {};
    const bIdx = CYC.map((c, i) => (c >= base0 && c <= base1) ? i : -1).filter(i => i >= 0);
    const sds = [], derivMeds = [];
    for (const w of wellsOf(run)) {
      if (!hasData(ch, w)) continue;
      const r = res[w] = { cq: null, tm: null, peaks: [], corr: null, deriv: null };
      const y = ch.amp[w];
      if (y && CYC.length) {
        const pts = bIdx.filter(i => y[i] != null);
        const f = linfit(pts.map(i => CYC[i]), pts.map(i => y[i]));
        r.fit = f;
        r.corr = y.map((v, i) => v == null ? null : v - (f.a + f.b * CYC[i]));
        const s = sd(pts.map(i => r.corr[i]));
        if (isFinite(s)) sds.push(s);
      }
      const m = ch.melt[w];
      if (m && TEMPS.length > 4) {
        r.deriv = negDeriv(m, TEMPS);
        const vals = r.deriv.filter(v => v != null).map(Math.abs);
        if (vals.length) derivMeds.push(median(vals));
      }
    }
    const autoThr = Math.max(10 * mean(sds), 1e-6) || 1;
    const thr = opts.thr ?? autoThr;
    const autoPeak = Math.max(10 * median(derivMeds), 1e-6) || 1;
    const pthr = opts.peakThr ?? autoPeak;
    for (const w in res) {
      const r = res[w];
      if (r.corr) r.cq = findCq(r.corr, thr, CYC);
      if (r.deriv) { r.peaks = findPeaks(r.deriv, pthr, TEMPS); r.tm = r.peaks.length ? r.peaks[0].t : null; }
    }
    return { res, thr, autoThr, pthr, autoPeak, base0, base1 };
  }

  // "77.1°C (12%), ..." for every melt peak after the main one
  function otherPeaks(peaks) {
    const p = peaks || [];
    if (p.length < 2) return "";
    return p.slice(1).map(q => `${q.t.toFixed(1)}°C (${Math.round(100 * q.h / p[0].h)}%)`).join(", ");
  }

  return { ROWS, wellsOf, mean, sd, median, hasData, defaultBaseline, linfit, findCq, negDeriv, findPeaks, analyze, otherPeaks };
})();
if (typeof module !== "undefined" && module.exports) module.exports = QPCR;
