import { test } from "node:test";
import assert from "node:assert/strict";
import "../js/share.js";
const { encode, decode, thin } = globalThis.QWaveShare;

const state = {
  scene: "doubleSlit", N: 256, k0: 4.2, sigma: 2.2, g: 0, display: "phase", speed: 8, detectorOn: true, detectorX: 1.5,
  packet: { center: [-11.25, 0.5], dir: [3, 1] },
  strokes: [{ erase: false, r: 1, h: 100, pts: [[1, 2], [3.05, -4.4]] }, { clear: true }, { erase: true, r: 0.5, h: 100, pts: [[0, 0]] }],
};

test("share links round-trip the app state", () => {
  const code = encode(state);
  assert.match(code, /^[\w-]+$/, "URL-safe characters only");
  const back = decode(code);
  assert.equal(back.scene, "doubleSlit"); assert.equal(back.N, 256); assert.equal(back.k0, 4.2);
  assert.equal(back.display, "phase"); assert.equal(back.detectorOn, true); assert.equal(back.detectorX, 1.5);
  assert.deepEqual(back.packet, state.packet);
  assert.deepEqual(back.strokes, state.strokes);
});

test("detector off and no custom packet", () => {
  const back = decode(encode({ ...state, detectorOn: false, packet: null, strokes: [] }));
  assert.equal(back.detectorOn, false); assert.equal(back.packet, null); assert.deepEqual(back.strokes, []);
});

test("malformed links are rejected", () => {
  for (const bad of ["", "abc", "e30", encode(state).slice(0, 20)]) assert.equal(decode(bad), null, bad);
});

test("thin keeps the ends and drops close points", () => {
  const pts = [[0, 0], [0.1, 0], [0.2, 0], [1, 0], [1.05, 0], [2, 0]];
  assert.deepEqual(thin(pts, 0.5), [[0, 0], [1, 0], [2, 0]]);
});
