const BactflowMeters = {
  history: { cpu: [], ram: [] },
  maxPoints: 40,
  timer: null,

  appRoot() {
    const el = document.getElementById("bf-meters");
    const root = (el && el.getAttribute("data-app-root")) || "";
    return String(root).replace(/\/$/, "");
  },

  apiUrl(path) {
    const root = this.appRoot();
    const clean = path.startsWith("/") ? path : `/${path}`;
    return `${root}${clean}`;
  },

  init() {
    if (!document.getElementById("bf-meters")) {
      return;
    }
    this.tick();
    this.timer = setInterval(() => this.tick(), 1000);
    document.addEventListener("visibilitychange", () => {
      if (!document.hidden) {
        this.tick();
      }
    });
  },

  level(pct) {
    if (pct >= 90) {
      return "is-hot";
    }
    if (pct >= 75) {
      return "is-warn";
    }
    return "";
  },

  push(kind, value) {
    const series = this.history[kind];
    series.push(Math.max(0, Math.min(100, Number(value) || 0)));
    if (series.length > this.maxPoints) {
      series.shift();
    }
  },

  drawSpark(svg, values, color) {
    if (!svg || values.length < 2) {
      return;
    }
    const w = 120;
    const h = 24;
    const max = 100;
    const pts = values.map((v, i) => {
      const x = (i / Math.max(values.length - 1, 1)) * w;
      const y = h - (v / max) * (h - 3) - 1.5;
      return `${x.toFixed(1)},${y.toFixed(1)}`;
    });
    const last = values[values.length - 1];
    const lx = ((values.length - 1) / Math.max(values.length - 1, 1)) * w;
    const ly = h - (last / max) * (h - 3) - 1.5;
    svg.innerHTML =
      `<polyline fill="none" stroke="${color}" stroke-width="1.6" points="${pts.join(" ")}"></polyline>` +
      `<circle cx="${lx.toFixed(1)}" cy="${ly.toFixed(1)}" r="2.2" fill="${color}"></circle>`;
  },

  setBar(fillEl, jobEl, pct, jobPct) {
    if (!fillEl) {
      return;
    }
    const host = Math.max(0, Math.min(100, Number(pct) || 0));
    fillEl.style.width = `${host}%`;
    fillEl.className = "bf-meter-fill " + this.level(host);
    if (jobEl) {
      const job = Math.max(0, Math.min(100, Number(jobPct) || 0));
      jobEl.style.width = `${job}%`;
    }
  },

  markWaiting(el, msg) {
    if (!el) {
      return;
    }
    el.classList.add("is-missing");
    el.textContent = msg || "waiting…";
  },

  apply(s) {
    if (!s || (s.ram_percent == null && s.cpu_percent == null)) {
      return false;
    }
    const cpuFill = document.getElementById("bf-cpu-fill");
    const ramFill = document.getElementById("bf-ram-fill");
    const cpuVal = document.getElementById("bf-cpu-val");
    const ramVal = document.getElementById("bf-ram-val");
    if (!cpuFill || !ramFill) {
      return false;
    }
    this.push("cpu", s.cpu_percent);
    this.push("ram", s.ram_percent);
    this.setBar(cpuFill, document.getElementById("bf-cpu-job"), s.cpu_percent, s.job_cpu_percent);
    this.setBar(
      ramFill,
      document.getElementById("bf-ram-job"),
      s.ram_percent,
      s.ram_total_gb ? (Number(s.job_rss_gb) / Number(s.ram_total_gb)) * 100 : 0
    );
    if (cpuVal) {
      cpuVal.classList.remove("is-missing");
      const jobBit = s.running && s.job_cores != null
        ? ` · <span class="job">job ${s.job_cores}/${s.cpu_count} cores</span>`
        : "";
      cpuVal.innerHTML = `${Number(s.cpu_percent || 0).toFixed(1)}%${jobBit}`;
    }
    if (ramVal) {
      ramVal.classList.remove("is-missing");
      const jobBit = s.running && s.job_rss_gb != null
        ? ` · <span class="job">job ${s.job_rss_gb} GB</span>`
        : "";
      ramVal.innerHTML = `${s.ram_used_gb} / ${s.ram_total_gb} GB (${Number(s.ram_percent || 0).toFixed(0)}%)${jobBit}`;
    }
    this.drawSpark(document.getElementById("bf-cpu-spark"), this.history.cpu, "#4fc1ff");
    this.drawSpark(document.getElementById("bf-ram-spark"), this.history.ram, "#c586c0");
    return true;
  },

  async tick() {
    if (document.hidden) {
      return;
    }
    const cpuVal = document.getElementById("bf-cpu-val");
    const ramVal = document.getElementById("bf-ram-val");
    try {
      let payload = null;
      const stamp = Date.now();
      const urls = [
        this.apiUrl("/sys_stats"),
        this.apiUrl("/bactflow_status"),
        this.apiUrl(`/static/sys_stats.json?t=${stamp}`),
        "sys_stats",
        `static/sys_stats.json?t=${stamp}`,
      ];
      for (const url of urls) {
        try {
          const res = await fetch(url, { cache: "no-store" });
          if (!res.ok) {
            continue;
          }
          const data = await res.json();
          if (data && (data.cpu_percent != null || data.ram_percent != null)) {
            payload = data;
            break;
          }
        } catch (err) {
          continue;
        }
      }
      if (!this.apply(payload)) {
        throw new Error("no stats fields");
      }
    } catch (err) {
      this.markWaiting(cpuVal);
      this.markWaiting(ramVal);
    }
  },
};

window.BactflowMeters = BactflowMeters;

document.addEventListener("DOMContentLoaded", () => {
  BactflowMeters.init();
});
