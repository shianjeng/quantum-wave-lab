/*
 * 3D surface view (WebGL2): the view grid as a lit height field, height = |ψ| (walls as low ridges),
 * colour = the same pixels as the 2D view (colour map, phase colours, wall tint).
 * It renders into its own offscreen canvas; the app copies that onto the visible canvas, so overlays,
 * PNG export and video recording work unchanged. Plain script: defines globalThis.QWaveSurface.
 */
(function (root) {
  "use strict";

  const VS = `#version 300 es
  precision highp float;
  in vec2 aPos;                                   // texel centres in [0, 1]²
  uniform sampler2D uHeight;
  uniform mat4 uMVP;
  uniform float uScale;
  uniform vec2 uTexel;
  out vec2 vUV;
  out vec3 vNormal;
  out float vH;
  void main() {
    float h = texture(uHeight, aPos).r;
    float hl = texture(uHeight, aPos - vec2(uTexel.x, 0.0)).r, hr = texture(uHeight, aPos + vec2(uTexel.x, 0.0)).r;
    float hd = texture(uHeight, aPos - vec2(0.0, uTexel.y)).r, hu = texture(uHeight, aPos + vec2(0.0, uTexel.y)).r;
    // world x, y in [-1, 1]: one texel is 2·texel wide, so dz/dx = (hr - hl)·scale / (4·texel)
    vNormal = normalize(vec3((hl - hr) * uScale / (4.0 * uTexel.x), (hd - hu) * uScale / (4.0 * uTexel.y), 1.0));
    vUV = aPos;
    vH = h;
    gl_Position = uMVP * vec4(aPos * 2.0 - 1.0, h * uScale, 1.0);
  }`;

  const FS = `#version 300 es
  precision highp float;
  in vec2 vUV;
  in vec3 vNormal;
  in float vH;
  uniform sampler2D uColor;                       // rows top to bottom, like ImageData
  uniform vec3 uLight;
  uniform vec3 uView;
  out vec4 outColor;
  void main() {
    vec3 base = texture(uColor, vec2(vUV.x, 1.0 - vUV.y)).rgb;
    vec3 n = normalize(vNormal), l = normalize(uLight), v = normalize(uView);
    float diff = max(dot(n, l), 0.0);
    float spec = pow(max(dot(n, normalize(l + v)), 0.0), 40.0);
    vec3 lit = base * (0.42 + 0.78 * diff) + vec3(1.0) * spec * (0.12 + 0.5 * vH);
    outColor = vec4(lit, 1.0);
  }`;

  // column-major 4×4 matrices
  function perspective(fovy, aspect, near, far) {
    const f = 1 / Math.tan(fovy / 2), nf = 1 / (near - far);
    return new Float32Array([f / aspect, 0, 0, 0, 0, f, 0, 0, 0, 0, (far + near) * nf, -1, 0, 0, 2 * far * near * nf, 0]);
  }
  function lookAt(eye, target, up) {
    const z = norm(sub(eye, target)), x = norm(cross(up, z)), y = cross(z, x);
    return new Float32Array([x[0], y[0], z[0], 0, x[1], y[1], z[1], 0, x[2], y[2], z[2], 0,
      -dot(x, eye), -dot(y, eye), -dot(z, eye), 1]);
  }
  function mul(a, b) {
    const o = new Float32Array(16);
    for (let c = 0; c < 4; c++) for (let r = 0; r < 4; r++) {
      let s = 0;
      for (let k = 0; k < 4; k++) s += a[k * 4 + r] * b[c * 4 + k];
      o[c * 4 + r] = s;
    }
    return o;
  }
  const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
  const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
  const norm = (a) => { const l = Math.hypot(...a) || 1; return [a[0] / l, a[1] / l, a[2] / l]; };

  class Surface3D {
    constructor(size = 800) {
      this.canvas = document.createElement("canvas");
      this.canvas.width = this.canvas.height = size;
      const gl = this.canvas.getContext("webgl2", { antialias: true, preserveDrawingBuffer: true });
      this.ok = Boolean(gl);
      if (!gl) return;
      this.gl = gl;
      this.prog = this._program(VS, FS);
      this.loc = {
        aPos: gl.getAttribLocation(this.prog, "aPos"),
        ...Object.fromEntries(["uHeight", "uColor", "uMVP", "uScale", "uTexel", "uLight", "uView"].map((u) => [u, gl.getUniformLocation(this.prog, u)])),
      };
      this.hTex = this._texture();
      this.cTex = this._texture();
      this.n = 0;
      this.camera = { yaw: -0.55, pitch: 0.78, dist: 3.1 };
    }

    _program(vs, fs) {
      const gl = this.gl, p = gl.createProgram();
      for (const [type, src] of [[gl.VERTEX_SHADER, vs], [gl.FRAGMENT_SHADER, fs]]) {
        const sh = gl.createShader(type);
        gl.shaderSource(sh, src); gl.compileShader(sh);
        if (!gl.getShaderParameter(sh, gl.COMPILE_STATUS)) throw new Error(gl.getShaderInfoLog(sh));
        gl.attachShader(p, sh);
      }
      gl.linkProgram(p);
      if (!gl.getProgramParameter(p, gl.LINK_STATUS)) throw new Error(gl.getProgramInfoLog(p));
      return p;
    }

    _texture() {
      const gl = this.gl, t = gl.createTexture();
      gl.bindTexture(gl.TEXTURE_2D, t);
      for (const [k, v] of [[gl.TEXTURE_MIN_FILTER, gl.LINEAR], [gl.TEXTURE_MAG_FILTER, gl.LINEAR],
        [gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE], [gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE]]) gl.texParameteri(gl.TEXTURE_2D, k, v);
      return t;
    }

    /** One vertex per texel centre, two triangles per cell. */
    _mesh(n) {
      const gl = this.gl, pos = new Float32Array(n * n * 2), idx = new Uint32Array((n - 1) * (n - 1) * 6);
      for (let j = 0; j < n; j++) for (let i = 0; i < n; i++) { pos[2 * (j * n + i)] = (i + 0.5) / n; pos[2 * (j * n + i) + 1] = (j + 0.5) / n; }
      let k = 0;
      for (let j = 0; j < n - 1; j++) for (let i = 0; i < n - 1; i++) {
        const a = j * n + i, b = a + 1, c = a + n, d = c + 1;
        idx[k++] = a; idx[k++] = b; idx[k++] = d; idx[k++] = a; idx[k++] = d; idx[k++] = c;
      }
      if (!this.vao) { this.vao = gl.createVertexArray(); this.vbo = gl.createBuffer(); this.ibo = gl.createBuffer(); }
      gl.bindVertexArray(this.vao);
      gl.bindBuffer(gl.ARRAY_BUFFER, this.vbo); gl.bufferData(gl.ARRAY_BUFFER, pos, gl.STATIC_DRAW);
      gl.enableVertexAttribArray(this.loc.aPos);
      gl.vertexAttribPointer(this.loc.aPos, 2, gl.FLOAT, false, 0, 0);
      gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, this.ibo); gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, idx, gl.STATIC_DRAW);
      gl.bindVertexArray(null);
      this.count = idx.length;
      this.n = n;
    }

    /** heights: Uint8Array n·n, row 0 = bottom (y = -20); rgba: n·n·4, row 0 = top (as ImageData). */
    draw(n, heights, rgba) {
      const gl = this.gl;
      if (n !== this.n) this._mesh(n);
      gl.pixelStorei(gl.UNPACK_ALIGNMENT, 1);
      gl.bindTexture(gl.TEXTURE_2D, this.hTex);
      gl.texImage2D(gl.TEXTURE_2D, 0, gl.R8, n, n, 0, gl.RED, gl.UNSIGNED_BYTE, heights);
      gl.bindTexture(gl.TEXTURE_2D, this.cTex);
      gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA8, n, n, 0, gl.RGBA, gl.UNSIGNED_BYTE, new Uint8Array(rgba.buffer, rgba.byteOffset, rgba.length));

      const { yaw, pitch, dist } = this.camera, target = [0, 0, 0.12];
      const eye = [target[0] + dist * Math.cos(pitch) * Math.sin(yaw), target[1] - dist * Math.cos(pitch) * Math.cos(yaw), target[2] + dist * Math.sin(pitch)];
      const mvp = mul(perspective((36 * Math.PI) / 180, 1, 0.05, 20), lookAt(eye, target, [0, 0, 1]));

      gl.viewport(0, 0, this.canvas.width, this.canvas.height);
      gl.clearColor(0.02, 0.027, 0.05, 1);
      gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
      gl.enable(gl.DEPTH_TEST);
      gl.useProgram(this.prog);
      gl.activeTexture(gl.TEXTURE0); gl.bindTexture(gl.TEXTURE_2D, this.hTex); gl.uniform1i(this.loc.uHeight, 0);
      gl.activeTexture(gl.TEXTURE1); gl.bindTexture(gl.TEXTURE_2D, this.cTex); gl.uniform1i(this.loc.uColor, 1);
      gl.uniformMatrix4fv(this.loc.uMVP, false, mvp);
      gl.uniform1f(this.loc.uScale, 0.55);
      gl.uniform2f(this.loc.uTexel, 1 / n, 1 / n);
      gl.uniform3f(this.loc.uLight, -0.4, -0.7, 1.0);
      gl.uniform3fv(this.loc.uView, sub(eye, target));
      gl.bindVertexArray(this.vao);
      gl.drawElements(gl.TRIANGLES, this.count, gl.UNSIGNED_INT, 0);
      gl.bindVertexArray(null);
      return this.canvas;
    }

    orbit(dx, dy) {
      const c = this.camera;
      c.yaw -= dx * 0.008;
      c.pitch = Math.min(1.45, Math.max(0.12, c.pitch + dy * 0.006));
    }
    zoom(f) { this.camera.dist = Math.min(7, Math.max(1.4, this.camera.dist * f)); }
  }

  root.QWaveSurface = { Surface3D };
})(globalThis);
