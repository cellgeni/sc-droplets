/* Interactive figures of the Explained pages: plain JavaScript and SVG, no libraries.
 *
 * A page places   <div class="scd-widget" data-widget="chord|cobweb|toy" data-pbmc="..." data-nuclei="..."></div>
 * and loads this file. The data files (docs/assets/explained/curves/*.json, written by
 * docs/_scripts/explained.py) hold each sample's totals, sorted and run-length coded.
 *
 * These widgets illustrate; they are not the methods. The chord and cobweb widgets compute the knee
 * from the identity the pages derive (on the cumulative curve, the farthest point from a chord is
 * the last barcode above the window's mean), which tests/test_explained.py checks against the
 * methods' own code. The page gives the package's own numbers beside them.
 */
(function () {
  "use strict";

  var SVGNS = "http://www.w3.org/2000/svg";
  var C = { blue: "#0072B2", vermillion: "#D55E00", green: "#009E73", purple: "#CC79A7", orange: "#E69F00", sky: "#56B4E9", grey: "#8c8c8c", light: "#d9d9d9", ink: "#27262b" };

  var STYLE = [
    ".scd-widget{border:1px solid #e6e1e8;border-radius:6px;padding:12px 14px;margin:1.2em 0;background:#fbfafc;font-size:14px;line-height:1.4}",
    ".scd-widget svg{display:block;width:100%;height:auto;overflow:visible}",
    ".scd-widget .scd-controls{display:flex;flex-wrap:wrap;gap:6px 16px;align-items:center;margin-bottom:8px}",
    ".scd-widget .scd-controls label{display:flex;flex-wrap:wrap;align-items:center;gap:2px 6px;max-width:100%}",
    ".scd-widget input[type=range]{width:140px;max-width:45vw}",
    ".scd-widget .scd-buttons{display:flex;flex-wrap:wrap;gap:4px}",
    ".scd-widget button{font:inherit;font-size:13px;padding:2px 8px;border:1px solid #cfc8d4;border-radius:4px;background:#fff;color:" + C.ink + ";cursor:pointer}",
    ".scd-widget button[aria-pressed=true]{background:#e9e4f2;border-color:#7253ed}",
    ".scd-widget .scd-readout{margin-top:6px;font-size:13px;color:#4a4a4a}",
    ".scd-widget .scd-readout b{color:" + C.ink + "}",
    ".scd-widget text{font-size:11px;fill:#4a4a4a}",
    ".scd-widget .scd-title{font-size:12px;fill:" + C.ink + "}",
  ].join("\n");

  function injectStyle() {
    if (document.getElementById("scd-widget-style")) return;
    var style = document.createElement("style");
    style.id = "scd-widget-style";
    style.textContent = STYLE;
    document.head.appendChild(style);
  }

  // -------------------------------------------------------------------------
  // Small helpers
  // -------------------------------------------------------------------------

  function svgEl(tag, attrs, parent) {
    var node = document.createElementNS(SVGNS, tag);
    for (var key in attrs) node.setAttribute(key, attrs[key]);
    if (parent) parent.appendChild(node);
    return node;
  }

  function htmlEl(tag, attrs, parent, text) {
    var node = document.createElement(tag);
    for (var key in attrs) node.setAttribute(key, attrs[key]);
    if (text !== undefined) node.textContent = text;
    if (parent) parent.appendChild(node);
    return node;
  }

  function fmt(x) {
    if (!isFinite(x)) return String(x);
    if (Math.abs(x) >= 100) return Math.round(x).toLocaleString("en-GB");
    if (Math.abs(x) >= 10) return x.toFixed(1);
    return x.toFixed(2);
  }

  function scale(d0, d1, r0, r1, log) {
    var f = log ? Math.log10 : function (v) { return v; };
    var a = f(d0), b = f(d1);
    var s = function (v) { return r0 + ((f(v) - a) / (b - a || 1)) * (r1 - r0); };
    s.d0 = d0; s.d1 = d1; s.r0 = r0; s.r1 = r1; s.log = log;
    return s;
  }

  function ticks(s) {
    var lo = Math.min(s.d0, s.d1), hi = Math.max(s.d0, s.d1), out = [];
    if (s.log) {
      for (var e = Math.ceil(Math.log10(lo) - 1e-9); e <= Math.floor(Math.log10(hi) + 1e-9); e++) out.push(Math.pow(10, e));
      return out;
    }
    var span = hi - lo, step = Math.pow(10, Math.floor(Math.log10(span / 4)));
    if (span / step > 10) step *= 5; else if (span / step > 5) step *= 2;
    for (var t = Math.ceil(lo / step) * step; t <= hi + step * 1e-9; t += step) out.push(t);
    return out;
  }

  function tickLabel(v, log) {
    if (log) {
      var e = Math.round(Math.log10(v));
      return e <= 4 ? Math.pow(10, e).toLocaleString("en-GB") : "1e" + e;
    }
    if (Math.abs(v) >= 1e6) return (v / 1e6).toLocaleString("en-GB") + "M";
    if (Math.abs(v) >= 1e4) return (v / 1e3).toLocaleString("en-GB") + "k";
    return v.toLocaleString("en-GB");
  }

  /** Axes with ticks and labels inside the box [x0, x1] x [y0, y1] (screen coordinates). */
  function drawAxes(g, sx, sy, xlabel, ylabel) {
    var x0 = Math.min(sx.r0, sx.r1), y1 = Math.max(sy.r0, sy.r1);
    svgEl("line", { x1: x0, x2: Math.max(sx.r0, sx.r1), y1: y1, y2: y1, stroke: "#999" }, g);
    svgEl("line", { x1: x0, x2: x0, y1: Math.min(sy.r0, sy.r1), y2: y1, stroke: "#999" }, g);
    ticks(sx).forEach(function (t) {
      var x = sx(t);
      svgEl("line", { x1: x, x2: x, y1: y1, y2: y1 + 4, stroke: "#999" }, g);
      svgEl("text", { x: x, y: y1 + 16, "text-anchor": "middle" }, g).textContent = tickLabel(t, sx.log);
    });
    ticks(sy).forEach(function (t) {
      var y = sy(t);
      svgEl("line", { x1: x0 - 4, x2: x0, y1: y, y2: y, stroke: "#999" }, g);
      svgEl("text", { x: x0 - 7, y: y + 4, "text-anchor": "end" }, g).textContent = tickLabel(t, sy.log);
    });
    svgEl("text", { x: (sx.r0 + sx.r1) / 2, y: y1 + 32, "text-anchor": "middle" }, g).textContent = xlabel;
    var yl = svgEl("text", { x: 0, y: 0, "text-anchor": "middle", transform: "translate(" + (x0 - 46) + "," + (sy.r0 + sy.r1) / 2 + ") rotate(-90)" }, g);
    yl.textContent = ylabel;
  }

  function polyline(points, attrs, parent) {
    var d = points.map(function (p) { return p[0].toFixed(1) + "," + p[1].toFixed(1); }).join(" ");
    var a = { points: d, fill: "none" };
    for (var key in attrs) a[key] = attrs[key];
    return svgEl("polyline", a, parent);
  }

  function clear(node) { while (node.firstChild) node.removeChild(node.firstChild); }

  function slider(parent, label, min, max, step, value, onInput) {
    var wrap = htmlEl("label", {}, parent);
    htmlEl("span", {}, wrap, label);
    var input = htmlEl("input", { type: "range", min: min, max: max, step: step, value: value }, wrap);
    var out = htmlEl("span", { "class": "scd-value" }, wrap);
    input.addEventListener("input", function () { onInput(+input.value); });
    return { input: input, out: out };
  }

  function buttons(parent, items, onPick) {
    var row = htmlEl("span", { "class": "scd-buttons", role: "group" }, parent);
    var list = items.map(function (item) {
      var b = htmlEl("button", { type: "button", "aria-pressed": "false" }, row, item.label);
      b.addEventListener("click", function () { onPick(item.value); });
      return { button: b, value: item.value };
    });
    return function select(value) {
      list.forEach(function (it) { it.button.setAttribute("aria-pressed", String(it.value === value)); });
    };
  }

  var cache = {};
  function load(url) {
    if (!cache[url]) cache[url] = fetch(url).then(function (r) { return r.json(); }).then(curve);
    return cache[url];
  }

  // -------------------------------------------------------------------------
  // A sample's sorted totals, from runs of equal counts
  // -------------------------------------------------------------------------

  function curve(data) {
    var nRuns = data.values.length, start = new Float64Array(nRuns + 1), before = new Float64Array(nRuns + 1);
    for (var j = 0; j < nRuns; j++) {
      start[j + 1] = start[j] + data.runs[j];
      before[j + 1] = before[j] + data.runs[j] * data.values[j];
    }
    var n = start[nRuns];
    function runOf(i) { // the run holding 0-based index i
      var lo = 0, hi = nRuns - 1;
      while (lo < hi) { var mid = (lo + hi + 1) >> 1; if (start[mid] <= i) lo = mid; else hi = mid - 1; }
      return lo;
    }
    function value(i) { return data.values[runOf(i)]; }
    function cum(i) { var j = runOf(i); return before[j] + (i - start[j] + 1) * data.values[j]; }
    /** The first 0-based index whose count is at most t (n if none). */
    function firstAtMost(t) {
      var lo = 0, hi = nRuns;
      while (lo < hi) { var mid = (lo + hi) >> 1; if (data.values[mid] <= t) hi = mid; else lo = mid + 1; }
      return lo < nRuns ? start[lo] : n;
    }
    /** UMI-tools' knee on the first m points, by the identity: barcodes 1..m-1 above their mean. */
    function knee(m) {
      if (m < 3) return 0;
      var mean = (cum(m - 1) - cum(0)) / (m - 1);
      return Math.max(0, Math.min(firstAtMost(mean), m) - 1);
    }
    /** Indices where the curve can turn: every run's first and last barcode, within [a, b]. */
    function corners(a, b) {
      var out = [a];
      for (var j = runOf(a); j < nRuns && start[j] <= b; j++) {
        if (start[j] > a) out.push(start[j]);
        var last = start[j + 1] - 1;
        if (last > a && last < b) out.push(last);
      }
      if (b > a) out.push(b);
      return out;
    }
    return { data: data, n: n, value: value, cum: cum, knee: knee, firstAtMost: firstAtMost, corners: corners };
  }

  // -------------------------------------------------------------------------
  // A: the chord
  // -------------------------------------------------------------------------

  function chordWidget(root) {
    var urls = { pbmc: root.dataset.pbmc, nuclei: root.dataset.nuclei };
    var state = { sample: "pbmc", view: "cumulative", a: 1, b: null, stretch: 1, curve: null };
    var controls = htmlEl("div", { "class": "scd-controls" }, root);
    var pickSample = buttons(controls, [{ label: "PBMC 1k", value: "pbmc" }, { label: "spinal cord nuclei", value: "nuclei" }], function (v) { state.sample = v; state.b = null; refresh(); });
    var pickView = buttons(controls, [{ label: "UMIs vs rank", value: "linear" }, { label: "log-log", value: "loglog" }, { label: "cumulative", value: "cumulative" }], function (v) { state.view = v; draw(); });
    var windows = htmlEl("div", { "class": "scd-controls" }, root);
    htmlEl("span", {}, windows, "The curve is drawn over:");
    var pickWindow = buttons(windows, [
      { label: "every barcode", value: "all" },
      { label: "the top 15,000 (dropbead)", value: "top15000" },
      { label: "≥ 1,000 UMIs (Kneedle)", value: "1000" },
      { label: "≥ 100 UMIs (zUMIs)", value: "100" },
    ], function (v) { setWindow(v); });
    var sliders = htmlEl("div", { "class": "scd-controls" }, root);
    var endA = slider(sliders, "first rank", 0, 1000, 1, 0, function (v) { state.a = toRank(v); if (state.a > state.b - 3) state.a = Math.max(1, state.b - 3); pickWindow(null); draw(); });
    var endB = slider(sliders, "last rank", 0, 1000, 1, 1000, function (v) { state.b = toRank(v); if (state.b < state.a + 3) state.b = state.a + 3; pickWindow(null); draw(); });
    var stretch = slider(sliders, "stretch the y axis", -1, 1, 0.01, 0, function (v) { state.stretch = Math.pow(2.5, v); draw(); });
    var svg = svgEl("svg", { viewBox: "0 0 640 360", role: "img", "aria-label": "A rank curve, a chord between two of its points, and the point farthest from the chord" }, root);
    var readout = htmlEl("div", { "class": "scd-readout", "aria-live": "polite" }, root);

    function toRank(v) { return Math.max(1, Math.round(Math.pow(10, (v / 1000) * Math.log10(state.curve.n)))); }
    function toSlider(r) { return Math.round((Math.log10(r) / Math.log10(state.curve.n)) * 1000); }

    function setWindow(kind) {
      var cv = state.curve;
      state.a = 1;
      if (kind === "all") state.b = cv.n;
      else if (kind === "top15000") state.b = Math.min(15000, cv.n);
      else state.b = Math.max(4, cv.firstAtMost(+kind - 1)); // barcodes with at least that many UMIs
      pickWindow(kind);
      draw();
    }

    function refresh() {
      pickSample(state.sample);
      load(urls[state.sample]).then(function (cv) {
        state.curve = cv;
        if (state.b === null) setWindow("all");
      });
    }

    function coords(i) {
      var cv = state.curve, r = i + 1, c = cv.value(i);
      if (state.view === "linear") return [r, c];
      if (state.view === "loglog") return [Math.log10(r), Math.log10(c)];
      return [r, cv.cum(i)];
    }

    function draw() {
      var cv = state.curve;
      if (!cv) return;
      pickView(state.view);
      endA.input.value = toSlider(state.a); endB.input.value = toSlider(state.b);
      endA.out.textContent = state.a.toLocaleString("en-GB"); endB.out.textContent = state.b.toLocaleString("en-GB");
      stretch.out.textContent = "×" + state.stretch.toFixed(2);
      var a = state.a - 1, b = state.b - 1;
      var idx = cv.corners(a, b), pts = idx.map(coords);
      var xs = pts.map(function (p) { return p[0]; }), ys = pts.map(function (p) { return p[1]; });
      var xmin = Math.min.apply(null, xs), xmax = Math.max.apply(null, xs), ymin = Math.min.apply(null, ys), ymax = Math.max.apply(null, ys);
      var W = 640, left = 70, right = 620, top = 16, h = Math.round(260 * state.stretch), H = top + h + 44;
      svg.setAttribute("viewBox", "0 0 " + W + " " + H);
      var sx = scale(xmin, xmax, left, right, false), sy = scale(ymin, ymax, top + h, top, false);
      clear(svg);
      var g = svgEl("g", {}, svg);
      var labels = { linear: ["rank", "UMIs"], loglog: ["log10 rank", "log10 UMIs"], cumulative: ["rank", "cumulative UMIs"] }[state.view];
      drawAxes(g, sx, sy, labels[0], labels[1]);
      var screen = pts.map(function (p) { return [sx(p[0]), sy(p[1])]; });
      var A = screen[0], B = screen[screen.length - 1];
      // the farthest point, measured on the screen: stretching cannot move it
      var best = 0, at = 0;
      for (var i = 0; i < screen.length; i++) {
        var area = Math.abs((B[0] - A[0]) * (screen[i][1] - A[1]) - (B[1] - A[1]) * (screen[i][0] - A[0]));
        if (area > best + 1e-9) { best = area; at = i; }
      }
      var P = screen[at];
      svgEl("polygon", { points: [A, B, P].map(function (p) { return p.join(","); }).join(" "), fill: C.sky, "fill-opacity": 0.18 }, g);
      polyline(screen, { stroke: C.grey, "stroke-width": 1.8 }, g);
      svgEl("line", { x1: A[0], y1: A[1], x2: B[0], y2: B[1], stroke: C.ink, "stroke-dasharray": "5 4" }, g);
      var ux = B[0] - A[0], uy = B[1] - A[1], len = Math.sqrt(ux * ux + uy * uy);
      ux /= len; uy /= len;
      var along = (P[0] - A[0]) * ux + (P[1] - A[1]) * uy, foot = [A[0] + along * ux, A[1] + along * uy];
      svgEl("line", { x1: P[0], y1: P[1], x2: foot[0], y2: foot[1], stroke: C.vermillion, "stroke-dasharray": "2 3" }, g);
      svgEl("line", { x1: P[0] - 90 * ux, y1: P[1] - 90 * uy, x2: P[0] + 90 * ux, y2: P[1] + 90 * uy, stroke: C.vermillion, "stroke-width": 1.5 }, g);
      svgEl("circle", { cx: P[0], cy: P[1], r: 5, fill: C.vermillion }, g);
      var rank = idx[at] + 1, text = "Farthest from the chord: <b>rank " + rank.toLocaleString("en-GB") + "</b> (" + cv.value(idx[at]).toLocaleString("en-GB") + " UMIs). 10x called " + cv.data.tenx.toLocaleString("en-GB") + ".";
      if (state.view === "cumulative") {
        var mean = (cv.cum(b) - cv.cum(a)) / (b - a), last = cv.firstAtMost(mean) - 1;
        text += " The chord's slope is the mean count of ranks " + (state.a + 1).toLocaleString("en-GB") + "–" + state.b.toLocaleString("en-GB") + ": <b>" + fmt(mean) + " UMIs</b>; the last barcode above it is rank " + (last + 1).toLocaleString("en-GB") + ".";
      }
      readout.innerHTML = text;
    }

    pickView(state.view);
    refresh();
  }

  // -------------------------------------------------------------------------
  // B: the refit loop as a map
  // -------------------------------------------------------------------------

  function cobwebWidget(root) {
    var urls = { pbmc: root.dataset.pbmc, nuclei: root.dataset.nuclei };
    var state = { sample: root.dataset.sample || "nuclei", multiplier: 3, start: null, curve: null };
    var controls = htmlEl("div", { "class": "scd-controls" }, root);
    var pickSample = buttons(controls, [{ label: "PBMC 1k", value: "pbmc" }, { label: "spinal cord nuclei", value: "nuclei" }], function (v) { state.sample = v; state.start = null; refresh(); });
    var mult = slider(controls, "window = k ×", 2, 8, 0.1, 3, function (v) { state.multiplier = v; draw(); });
    var starts = htmlEl("div", { "class": "scd-controls" }, root);
    var first = slider(starts, "first knee", 0, 1000, 1, 1000, function (v) { state.start = Math.max(2, Math.round(Math.pow(10, 1 + (v / 1000) * (Math.log10(state.curve.n) - 1)))); draw(); });
    var reset = htmlEl("button", { type: "button" }, starts, "from every barcode");
    reset.addEventListener("click", function () { state.start = null; draw(); });
    var svg = svgEl("svg", { viewBox: "0 0 640 400", role: "img", "aria-label": "The map from one knee to the next, the diagonal, and the passes drawn as a cobweb" }, root);
    var readout = htmlEl("div", { "class": "scd-readout", "aria-live": "polite" }, root);

    function refresh() {
      pickSample(state.sample);
      load(urls[state.sample]).then(function (cv) { state.curve = cv; draw(); });
    }

    function draw() {
      var cv = state.curve;
      if (!cv) return;
      var c = state.multiplier, n = cv.n;
      mult.out.textContent = c.toFixed(1);
      function next(k) { return cv.knee(Math.min(Math.floor(c * k), n)); }
      var k0 = state.start === null ? cv.knee(n) : state.start;
      first.input.value = Math.round(((Math.log10(Math.max(10, k0)) - 1) / (Math.log10(n) - 1)) * 1000);
      first.out.textContent = k0.toLocaleString("en-GB");
      var lo = 10, hi = n;
      var W = 640, left = 70, right = 620, top = 12, bottom = 350;
      var sx = scale(lo, hi, left, right, true), sy = scale(lo, hi, bottom, top, true);
      clear(svg);
      var g = svgEl("g", {}, svg);
      drawAxes(g, sx, sy, "knee k (log)", "next knee (log)");
      var grid = [], map = [];
      for (var i = 0; i <= 400; i++) {
        var k = Math.round(Math.pow(10, Math.log10(lo) + (i / 400) * (Math.log10(hi) - Math.log10(lo))));
        if (grid.length && k === grid[grid.length - 1]) continue;
        grid.push(k); map.push(next(k));
      }
      svgEl("line", { x1: sx(lo), y1: sy(lo), x2: sx(hi), y2: sy(hi), stroke: C.ink, "stroke-dasharray": "5 4", "stroke-width": 0.8 }, g);
      polyline(grid.map(function (k, j) { return [sx(k), sy(Math.max(lo, map[j]))]; }), { stroke: C.grey, "stroke-width": 2 }, g);
      var tenx = cv.data.tenx;
      svgEl("line", { x1: sx(tenx), x2: sx(tenx), y1: top, y2: bottom, stroke: C.blue, "stroke-dasharray": "2 3" }, g);
      svgEl("text", { x: sx(tenx) + 4, y: top + 10, fill: C.blue }, g).textContent = "10x";
      // fixed points on the grid
      var fixed = [];
      for (var j = 1; j < grid.length; j++) {
        var before = Math.sign(map[j - 1] - grid[j - 1]), after = Math.sign(map[j] - grid[j]);
        if (before !== after && before !== 0 && grid[j] >= 500) fixed.push({ k: grid[j], stable: before > 0 });
      }
      fixed.forEach(function (f) {
        svgEl("circle", { cx: sx(f.k), cy: sy(f.k), r: 5.5, fill: f.stable ? C.ink : "#fff", stroke: C.ink, "stroke-width": 1.5 }, g);
      });
      // the passes
      var path = [[k0, k0]], k = k0, passes = 1, collapsed = false;
      for (var p = 0; p < 200; p++) {
        var nk = next(k);
        if (nk < 1) { collapsed = true; break; }
        path.push([k, nk], [nk, nk]);
        passes++;
        if (nk === k) break;
        k = nk;
      }
      polyline(path.map(function (q) { return [sx(Math.max(lo, q[0])), sy(Math.max(lo, q[1]))]; }), { stroke: C.vermillion, "stroke-width": 1.5 }, g);
      svgEl("circle", { cx: sx(Math.max(lo, k0)), cy: sy(Math.max(lo, k0)), r: 4, fill: "#fff", stroke: C.vermillion, "stroke-width": 1.5 }, g);
      var stable = fixed.filter(function (f) { return f.stable; }).length;
      var text = collapsed
        ? "With a window of k × " + c.toFixed(1) + " the knee shrinks pass after pass until no barcode is left: <b>the loop collapses</b>."
        : "After " + passes + " passes the loop settles at <b>" + (k + 1).toLocaleString("en-GB") + " cells</b>.";
      text += " The map crosses the diagonal at " + fixed.length + " point" + (fixed.length === 1 ? "" : "s") + " (" + stable + " stable, filled). 10x called " + tenx.toLocaleString("en-GB") + ".";
      readout.innerHTML = text;
    }

    refresh();
  }

  // -------------------------------------------------------------------------
  // C: a toy rank curve
  // -------------------------------------------------------------------------

  function mulberry32(seed) {
    return function () {
      seed |= 0; seed = (seed + 0x6D2B79F5) | 0;
      var t = Math.imul(seed ^ (seed >>> 15), 1 | seed);
      t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
  }

  function lognormals(n, median, sigma, seed) {
    var rand = mulberry32(seed), out = new Float64Array(n), mu = Math.log(median);
    for (var i = 0; i < n; i++) {
      var u = 1 - rand(), v = rand();
      var z = Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v);
      out[i] = Math.max(1, Math.round(Math.exp(mu + sigma * z)));
    }
    return out;
  }

  function toyWidget(root) {
    var state = { cells: 1500, cellSize: 3000, empty: 8000, emptySize: 40 };
    var controls = htmlEl("div", { "class": "scd-controls" }, root);
    function logSlider(label, key, min, max) {
      var s = slider(controls, label, Math.log10(min), Math.log10(max), 0.01, Math.log10(state[key]), function (v) {
        state[key] = Math.round(Math.pow(10, v)); schedule();
      });
      s.key = key;
      return s;
    }
    var all = [logSlider("cells", "cells", 100, 20000), logSlider("UMIs per cell (median)", "cellSize", 200, 30000),
               logSlider("empty droplets", "empty", 1000, 200000), logSlider("UMIs per empty droplet", "emptySize", 2, 1000)];
    var svg = svgEl("svg", { viewBox: "0 0 640 300", role: "img", "aria-label": "A simulated rank curve and the histogram of its log totals" }, root);
    var readout = htmlEl("div", { "class": "scd-readout", "aria-live": "polite" }, root);
    var pending = false;
    function schedule() { if (!pending) { pending = true; requestAnimationFrame(function () { pending = false; draw(); }); } }

    function draw() {
      all.forEach(function (s) { s.out.textContent = state[s.key].toLocaleString("en-GB"); });
      var groups = [
        { kind: 0, values: lognormals(state.cells, state.cellSize, 0.6, 1), colour: C.blue, label: "cells" },
        { kind: 1, values: lognormals(state.empty, state.emptySize, 0.5, 2), colour: C.orange, label: "empty droplets" },
        { kind: 2, values: lognormals(30000, 2, 0.8, 3), colour: C.grey, label: "error barcodes" },
      ];
      var total = groups.reduce(function (s, g) { return s + g.values.length; }, 0);
      var keys = new Float64Array(total), kinds = new Uint8Array(total), at = 0;
      groups.forEach(function (grp) { for (var i = 0; i < grp.values.length; i++) { keys[at] = grp.values[i]; kinds[at] = grp.kind; at++; } });
      // sort by count, decreasing, carrying the kind along
      var order = new Uint32Array(total);
      for (var i = 0; i < total; i++) order[i] = i;
      order.sort(function (x, y) { return keys[y] - keys[x]; });
      clear(svg);
      var g = svgEl("g", {}, svg);
      var maxCount = keys[order[0]];
      var sx = scale(1, total, 60, 330, true), sy = scale(1, maxCount * 1.2, 250, 14, true);
      drawAxes(g, sx, sy, "rank (log)", "UMIs (log)");
      var picks = {};
      for (var p = 0; p <= 1600; p++) picks[Math.min(total - 1, Math.round(Math.pow(10, (p / 1600) * Math.log10(total))) - 1)] = true;
      Object.keys(picks).forEach(function (key) {
        var r = +key, i = order[r];
        svgEl("circle", { cx: sx(r + 1), cy: sy(keys[i]), r: 1.6, fill: groups[kinds[i]].colour }, g);
      });
      // the histogram of log10 totals, stacked
      var bins = 50, top = Math.log10(maxCount) + 0.05, counts = groups.map(function () { return new Float64Array(bins); });
      for (var j = 0; j < total; j++) {
        var b = Math.min(bins - 1, Math.floor((Math.log10(keys[j]) / top) * bins));
        counts[kinds[j]][b]++;
      }
      var hmax = 1;
      for (var bb = 0; bb < bins; bb++) hmax = Math.max(hmax, counts[0][bb] + counts[1][bb] + counts[2][bb]);
      var hx = scale(0, top, 410, 620, false), hy = scale(1, hmax * 1.3, 250, 14, true);
      drawAxes(g, hx, hy, "UMIs (log10)", "barcodes (log)");
      var width = (hx(top) - hx(0)) / bins;
      for (var k = 0; k < bins; k++) {
        var base = 0;
        [2, 1, 0].forEach(function (kind) {
          var v = counts[kind][k];
          if (!v) return;
          var y0 = hy(Math.max(1, base)), y1 = hy(base + v);
          svgEl("rect", { x: hx((k / bins) * top), y: y1, width: Math.max(0.5, width - 0.6), height: Math.max(0.5, (base ? y0 : 250) - y1), fill: groups[kind].colour }, g);
          base += v;
        });
      }
      groups.forEach(function (grp, i) {
        svgEl("rect", { x: 330 - 110, y: 20 + i * 15, width: 9, height: 9, fill: grp.colour }, g);
        svgEl("text", { x: 330 - 97, y: 28 + i * 15 }, g).textContent = grp.label;
      });
      var ratio = state.cellSize / state.emptySize;
      readout.innerHTML = "Cells hold <b>" + (ratio >= 10 ? Math.round(ratio).toLocaleString("en-GB") : ratio.toFixed(1)) + "×</b> the UMIs of an empty droplet (medians). " +
        (ratio > 30 ? "The cliff is sharp and the two humps stand apart." : ratio > 8 ? "The cliff softens into a bend; the humps begin to touch." : "There is no cliff to find: the humps have merged.");
    }

    draw();
  }

  // -------------------------------------------------------------------------

  var WIDGETS = { chord: chordWidget, cobweb: cobwebWidget, toy: toyWidget };

  function start() {
    injectStyle();
    var nodes = document.querySelectorAll(".scd-widget[data-widget]");
    for (var i = 0; i < nodes.length; i++) {
      var make = WIDGETS[nodes[i].dataset.widget];
      if (make && !nodes[i].dataset.ready) {
        nodes[i].dataset.ready = "1";
        make(nodes[i]);
      }
    }
  }

  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", start);
  else start();
})();
