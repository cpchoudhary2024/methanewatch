(() => {
  "use strict";
  const svg = document.getElementById("study-map");
  const years = document.getElementById("years");
  const ns = "http://www.w3.org/2000/svg";
  let season = 0;
  let mode = "uniform";

  const el = (name, attrs = {}, parent = svg) => {
    const node = document.createElementNS(ns, name);
    for (const [key, value] of Object.entries(attrs)) node.setAttribute(key, String(value));
    parent.appendChild(node);
    return node;
  };
  const byId = id => document.getElementById(id);
  const pct = value => Math.round(value * 100) + "%";
  const est = value => value == null ? "n/a" : value.toFixed(2) + "×";
  const pText = value => value == null ? "p unavailable" : "p = " + value.toFixed(value < 0.01 ? 3 : 2);
  const xy = (point, bbox) => [
    32 + (point[1] - bbox[0]) / (bbox[2] - bbox[0]) * 656,
    450 - (point[0] - bbox[1]) / (bbox[3] - bbox[1]) * 416
  ];
  const text = (x, y, value, anchor = "start") => {
    const node = el("text", {x, y, "text-anchor": anchor, fill: "#6e8582", "font-size": 11, "font-family": "Arial,sans-serif"});
    node.textContent = value;
  };
  const dots = (points, bbox, style) => {
    const layer = el("g", style.group || {});
    for (const point of points) {
      const [cx, cy] = xy(point, bbox);
      el("circle", {cx, cy, r: style.radius, fill: style.fill, opacity: style.opacity}, layer);
    }
  };

  function drawMap(basin) {
    svg.replaceChildren();
    el("rect", {x: 0, y: 0, width: 720, height: 485, fill: "#f6f8f3"});
    for (let i = 0; i < 5; i++) {
      const x = 32 + i * 164, y = 34 + i * 104;
      el("line", {x1: x, y1: 34, x2: x, y2: 450, stroke: "#dce7dd", "stroke-width": 1});
      el("line", {x1: 32, y1: y, x2: 688, y2: y, stroke: "#dce7dd", "stroke-width": 1});
    }
    text(33, 472, "105° W");
    text(688, 472, "101° W", "end");
    text(15, 42, "33.5° N");
    text(15, 445, "30.8° N");
    dots(basin.flares, basin.bbox, {radius: 1.7, fill: "#d8833d", opacity: .48});
    dots(basin.controls[mode], basin.bbox, {radius: 2.7, fill: "#668998", opacity: .67});
    for (const candidate of basin.candidates) {
      const [cx, cy] = xy(candidate, basin.bbox);
      const circle = el("circle", {
        cx, cy, r: 6, fill: "#b74e40", stroke: "#fff", "stroke-width": 2,
        tabindex: 0, role: "button", "aria-label": candidate[3] + ", about " + candidate[2] + " tonnes per hour"
      });
      const label = candidate[3] + " · approximately " + candidate[2] + " tonnes CH₄ per hour";
      el("title", {}, circle).textContent = label;
      const select = () => { byId("selected-candidate").textContent = label + ". Candidate centre only; source unconfirmed."; };
      circle.addEventListener("click", select);
      circle.addEventListener("keydown", event => {
        if (event.key === "Enter" || event.key === " ") { event.preventDefault(); select(); }
      });
    }
  }

  function render() {
    const basin = DATA[season];
    if (!basin) return;
    years.querySelectorAll("button").forEach((button, i) => button.setAttribute("aria-pressed", String(i === season)));
    document.querySelectorAll("[data-mode]").forEach(button => button.setAttribute("aria-pressed", String(button.dataset.mode === mode)));
    byId("candidate-share").textContent = basin.n_candidates + " candidates / " + basin.n_overpasses + " overpasses";
    byId("candidate-percent").textContent = pct(basin.frac_det);
    byId("candidate-bar").style.width = pct(basin.frac_det);
    for (const key of ["uniform", "wells"]) {
      const item = basin.flare[key];
      byId(key + "-percent").textContent = pct(item.frac);
      byId(key + "-bar").style.width = pct(item.frac);
      byId(key + "-estimate").textContent = est(item.e);
      byId(key + "-detail").textContent = (item.ci ? "95% day interval " + item.ci[0].toFixed(2) + " to " + item.ci[1].toFixed(2) + "; " : "") + pText(item.p);
      const matched = basin.opportunity && basin.opportunity[key];
      byId("matched-" + key).textContent = matched && matched.enrichment != null
        ? est(matched.enrichment) + " (" + matched.ci_days[0].toFixed(2) + " to " + matched.ci_days[1].toFixed(2) + ")"
        : "not calculated";
    }
    byId("interpretation").textContent = "In " + basin.year + ", the box-wide uniform control gives enrichment. The well-density and QA-pixel matched day intervals include one. Candidate locations are identical in each comparison.";
    byId("selected-candidate").textContent = "Select a candidate to inspect its date and estimated rate.";
    drawMap(basin);
  }

  DATA.forEach((basin, i) => {
    const button = document.createElement("button");
    button.type = "button";
    button.textContent = basin.year;
    button.addEventListener("click", () => { season = i; render(); });
    years.appendChild(button);
  });
  document.querySelectorAll("[data-mode]").forEach(button => {
    button.addEventListener("click", () => { mode = button.dataset.mode; render(); });
  });
  render();
})();
