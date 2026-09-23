/*
 * Share links: the app state (scene, parameters, custom packet, drawn wall strokes) as a compact
 * base64url string for the URL hash. Plain script (globalThis.QWaveShare); pure, tested in Node.
 *
 * Coordinates are stored as integers in units of 0.05, which is finer than the finest grid (0.1).
 */
(function (root) {
  "use strict";
  const VERSION = 1, Q = 20;

  const q = (v) => Math.round(v * Q);
  const r2 = (v) => Math.round(v * 100) / 100;

  function toBase64Url(text) {
    const bytes = new TextEncoder().encode(text);
    let bin = "";
    for (const b of bytes) bin += String.fromCharCode(b);
    return btoa(bin).replace(/\+/g, "-").replace(/\//g, "_").replace(/=+$/, "");
  }

  function fromBase64Url(s) {
    const b64 = s.replace(/-/g, "+").replace(/_/g, "/") + "===".slice((s.length + 3) % 4);
    const bin = atob(b64);
    return new TextDecoder().decode(Uint8Array.from(bin, (c) => c.charCodeAt(0)));
  }

  /**
   * s: {scene, N, k0, sigma, g, display, speed, detectorOn, detectorX, packet: {center, dir} | null,
   *     strokes: [{clear: true} | {erase, r, h, pts: [[x, y], ...]}]}
   */
  function encode(s) {
    const o = {
      v: VERSION, sc: s.scene, n: s.N, k: r2(s.k0), s: r2(s.sigma), g: s.g, d: s.display, sp: s.speed,
      det: s.detectorOn ? r2(s.detectorX) : null,
      pk: s.packet ? [...s.packet.center.map(q), ...s.packet.dir.map(q)] : null,
      st: s.strokes.map((st) => (st.clear ? 0 : [st.erase ? 1 : 0, r2(st.r), Math.round(st.h), st.pts.flatMap(([x, y]) => [q(x), q(y)])])),
    };
    return toBase64Url(JSON.stringify(o));
  }

  const num = (v) => typeof v === "number" && Number.isFinite(v);

  /** Inverse of encode; returns null for anything malformed (the app then just ignores the link). */
  function decode(str) {
    try {
      const o = JSON.parse(fromBase64Url(str));
      if (o.v !== VERSION || typeof o.sc !== "string" || ![o.n, o.k, o.s, o.g, o.sp].every(num)) return null;
      const strokes = [];
      for (const st of o.st ?? []) {
        if (st === 0) { strokes.push({ clear: true }); continue; }
        const [erase, r, h, flat] = st;
        if (![r, h].every(num) || !Array.isArray(flat) || flat.length % 2 || !flat.every(num)) return null;
        const pts = [];
        for (let i = 0; i < flat.length; i += 2) pts.push([flat[i] / Q, flat[i + 1] / Q]);
        strokes.push({ erase: erase === 1, r, h, pts });
      }
      let packet = null;
      if (Array.isArray(o.pk)) {
        if (o.pk.length !== 4 || !o.pk.every(num)) return null;
        packet = { center: [o.pk[0] / Q, o.pk[1] / Q], dir: [o.pk[2] / Q, o.pk[3] / Q] };
      }
      return {
        scene: o.sc, N: o.n, k0: o.k, sigma: o.s, g: o.g, display: typeof o.d === "string" ? o.d : "amplitude",
        speed: o.sp, detectorOn: num(o.det), detectorX: num(o.det) ? o.det : 0, packet, strokes,
      };
    } catch (e) {
      return null;
    }
  }

  /** Drop points closer than minDist to the previous kept one (the last point is always kept). */
  function thin(pts, minDist) {
    if (pts.length <= 2) return pts.slice();
    const out = [pts[0]];
    for (let i = 1; i < pts.length - 1; i++) {
      const [x, y] = pts[i], [px, py] = out[out.length - 1];
      if (Math.hypot(x - px, y - py) >= minDist) out.push(pts[i]);
    }
    out.push(pts[pts.length - 1]);
    return out;
  }

  root.QWaveShare = { encode, decode, thin };
})(globalThis);
