const BactflowTerminal = {
  quiet: true,
  logEl: null,
  statusEl: null,
  maxLines: 1800,

  init() {
    this.logEl = document.getElementById("output-bactflow");
    this.statusEl = document.getElementById("bf-status");
    const toggle = document.getElementById("logVerbose");
    if (toggle) {
      this.quiet = !toggle.checked;
      toggle.addEventListener("change", () => {
        this.quiet = !toggle.checked;
      });
    }
  },

  stripAnsi(text) {
    return String(text || "")
      .replace(/\u001b\[[0-9;?]*[A-Za-z]/g, "")
      .replace(/\u001b\][^\u0007]*\u0007/g, "")
      .replace(/[\x00-\x08\x0b\x0c\x0e-\x1f]/g, "");
  },

  classify(line) {
    if (/ERROR|Error executing|terminated with an error|ModuleNotFoundError|Traceback|exit status|Command error/i.test(line)) {
      return "log-error";
    }
    if (/WARN|WARNING/i.test(line)) {
      return "log-warn";
    }
    if (/✔|successfully|Process completed|Cached process|fasta files are ready|started :\)/i.test(line)) {
      return "log-ok";
    }
    if (/\[100%\]/.test(line)) {
      return "log-ok";
    }
    if (/Launching|N E X T F L O W|Using Java|executor >/i.test(line)) {
      return "log-info";
    }
    if (/Submitted process|Cached process/.test(line) || (/\]\s+\S+/.test(line) && /\[\w{2}\/\w+\]/.test(line))) {
      return "log-process";
    }
    if (/^\[-\s*\]/.test(line) || /Tip: you can|Check '\.nextflow\.log'/.test(line)) {
      return "log-muted";
    }
    return "log-default";
  },

  isInteresting(line) {
    const t = line.trim();
    if (!t) {
      return false;
    }
    if (this.quiet === false) {
      return true;
    }
    if (/^\[-\s*\]/.test(t) && !/\[\d/.test(t)) {
      return false;
    }
    if (/Tip: you can|Check '\.nextflow\.log'|WORKFLOW OUTPUT DEFINITION|is available - Please consider/.test(t)) {
      return false;
    }
    return /ERROR|WARN|✔|Launching|Process completed|started|failed|ready in|Using Java|N E X T F L O W|executor >|Submitted process|Cached process|\[100%\]|\[  0%\]|Caused by:|Command exit status|running SPAdes|running Unicycler|running Flye|Running circlator|Stream disconnected|taxonomy|bakta|checkm|quast/i.test(t);
  },

  setStatus(text, kind) {
    if (!this.statusEl) {
      return;
    }
    this.statusEl.textContent = text;
    this.statusEl.className = "bf-status" + (kind ? " is-" + kind : "");
  },

  clear() {
    if (this.logEl) {
      this.logEl.replaceChildren();
    }
    this.setStatus("Ready", "");
  },

  append(raw, force) {
    if (!this.logEl) {
      this.init();
    }
    const line = this.stripAnsi(raw);
    if (!force && !this.isInteresting(line)) {
      return;
    }
    const row = document.createElement("div");
    row.className = "log-line " + this.classify(line);
    row.textContent = line || " ";
    this.logEl.appendChild(row);

    while (this.logEl.childNodes.length > this.maxLines) {
      this.logEl.removeChild(this.logEl.firstChild);
    }
    this.logEl.scrollTop = this.logEl.scrollHeight;

    if (/ERROR|terminated with an error|ModuleNotFoundError/i.test(line)) {
      this.setStatus(line.replace(/\s+/g, " ").slice(0, 140), "error");
    } else if (/\[100%\].*✔|Process completed|successfully/i.test(line)) {
      this.setStatus(line.replace(/\s+/g, " ").slice(0, 140), "ok");
    } else if (/Submitted process|running |Launching/i.test(line)) {
      this.setStatus(line.replace(/\s+/g, " ").slice(0, 140), "run");
    } else if (/WARN/i.test(line)) {
      this.setStatus(line.replace(/\s+/g, " ").slice(0, 140), "warn");
    }
  },

  appendMany(lines) {
    lines.forEach((line) => this.append(line));
  }
};

const BactflowMeters = {
  history: { cpu: [], ram: [] },
  maxPoints: 40,
  timer: null,

  init() {
    if (!document.getElementById("bf-meters")) {
      return;
    }
    this.tick();
    this.timer = setInterval(() => this.tick(), 1000);
    document.addEventListener("visibilitychange", () => {
      if (document.hidden) {
        return;
      }
      this.tick();
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
      `<polyline fill="none" stroke="${color}" stroke-width="1.6" stroke-linejoin="round" points="${pts.join(" ")}"></polyline>` +
      `<circle cx="${lx.toFixed(1)}" cy="${ly.toFixed(1)}" r="1.8" fill="${color}"></circle>`;
  },

  setBar(fillEl, jobEl, percent, jobPercent) {
    const pct = Math.max(0, Math.min(100, Number(percent) || 0));
    fillEl.style.width = pct.toFixed(1) + "%";
    fillEl.className = "bf-meter-fill " + this.level(pct);
    if (jobEl) {
      const job = Math.max(0, Math.min(100, Number(jobPercent) || 0));
      jobEl.style.width = job.toFixed(1) + "%";
    }
  },

  async tick() {
    if (document.hidden) {
      return;
    }
    const cpuFill = document.getElementById("bf-cpu-fill");
    const ramFill = document.getElementById("bf-ram-fill");
    const cpuVal = document.getElementById("bf-cpu-val");
    const ramVal = document.getElementById("bf-ram-val");
    if (!cpuFill || !ramFill) {
      return;
    }
    try {
      const res = await fetch("/sys_stats", { cache: "no-store" });
      if (!res.ok) {
        throw new Error("HTTP " + res.status);
      }
      const s = await res.json();
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
    } catch (err) {
      if (cpuVal) {
        cpuVal.classList.add("is-missing");
        cpuVal.textContent = "restart the app";
      }
      if (ramVal) {
        ramVal.classList.add("is-missing");
        ramVal.textContent = "no /sys_stats";
      }
    }
  }
};

document.addEventListener("DOMContentLoaded", () => {
  BactflowTerminal.init();
  BactflowMeters.init();
});
