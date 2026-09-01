const BactflowProcessEta = {
  processes: new Map(),
  hashToProc: new Map(),
  hashToName: new Map(),
  durations: new Map(),
  workflow: { startedAt: null, percent: 0, completed: 0, total: 0 },
  tickTimer: null,

  hintsSec: {
    envsetup: 240,
    testify: 45,
    coverage_filt: 180,
    readfilt: 180,
    concat_reads: 120,
    assembly_flye1: 1500,
    assembly_flye2: 1500,
    assembly_pacbio: 1500,
    assembly_spades: 2100,
    assembly_unicycler: 2700,
    circulator: 300,
    quast: 240,
    quastcheck: 240,
    bakta: 900,
    gtdbtk: 600,
    checkm: 600
  },

  reset() {
    this.processes.clear();
    this.hashToProc.clear();
    this.hashToName.clear();
    this.durations.clear();
    this.workflow = { startedAt: null, percent: 0, completed: 0, total: 0 };
    if (this.tickTimer) {
      clearInterval(this.tickTimer);
      this.tickTimer = null;
    }
  },

  ensureTick(terminal) {
    if (this.tickTimer) {
      return;
    }
    this.tickTimer = setInterval(() => this.refreshAll(terminal), 1000);
  },

  normalizeName(name) {
    return String(name || "").toLowerCase().replace(/[^a-z0-9]+/g, "");
  },

  hintFor(name) {
    const key = this.normalizeName(name);
    if (this.durations.has(key)) {
      return this.durations.get(key);
    }
    if (this.hintsSec[key]) {
      return this.hintsSec[key];
    }
    for (const [hintKey, sec] of Object.entries(this.hintsSec)) {
      if (key.includes(hintKey) || hintKey.includes(key)) {
        return sec;
      }
    }
    return 180;
  },

  stopTick() {
    if (this.tickTimer) {
      clearInterval(this.tickTimer);
      this.tickTimer = null;
    }
  },

  markDone(proc) {
    if (!proc || proc.done) {
      return;
    }
    proc.done = true;
    proc.finishedAt = Date.now();
    proc.frozenElapsed = (proc.finishedAt - proc.startedAt) / 1000;
    proc.percent = 100;
    this.durations.set(this.normalizeName(proc.name), proc.frozenElapsed);
  },

  markHashDone(hash, terminal) {
    const proc = this.hashToProc.get(hash);
    if (proc) {
      this.markDone(proc);
      if (proc.row && terminal) {
        this.renderRow(proc, terminal.classify(proc.lastLine));
      }
      return proc;
    }
    return null;
  },

  scanLineForCompletion(line, terminal) {
    const t = line.trim();
    if (!t || /\[Task monitor\]|TaskPollingMonitor|TaskHandler\[/i.test(t)) {
      return false;
    }

    let matched = false;
    const hashMatches = t.matchAll(/\[([0-9a-f/]+)\]/gi);
    for (const hm of hashMatches) {
      const hash = hm[1];
      if (
        /✔|COMPLETED|Task completed/i.test(t) ||
        /\[\s*100%\]/i.test(t) ||
        (/process\s+>/i.test(t) && /\[\s*\d+%\]/i.test(t) && /\b1\s+of\s+1\b/i.test(t))
      ) {
        if (this.markHashDone(hash, terminal)) {
          matched = true;
        }
      }
    }
    return matched;
  },

  finalizeAll(terminal, reason) {
    const now = Date.now();
    this.processes.forEach((proc) => {
      if (!proc.done) {
        proc.done = true;
        proc.finishedAt = now;
        proc.frozenElapsed = (now - proc.startedAt) / 1000;
        if (proc.row) {
          this.renderRow(proc, terminal.classify(proc.lastLine));
        }
      }
    });
    this.stopTick();
    if (reason === "completed") {
      terminal.setStatus("All processes completed", "ok");
    }
  },

  allDone() {
    return this.processes.size > 0 && [...this.processes.values()].every((p) => p.done);
  },

  formatDuration(seconds) {
    if (!Number.isFinite(seconds) || seconds < 0) {
      return "…";
    }
    if (seconds < 45) {
      return `${Math.max(1, Math.round(seconds))}s`;
    }
    if (seconds < 3600) {
      const m = Math.floor(seconds / 60);
      const s = Math.round(seconds % 60);
      return `${m}m ${s}s`;
    }
    const h = Math.floor(seconds / 3600);
    const m = Math.floor((seconds % 3600) / 60);
    return `${h}h ${m}m`;
  },

  formatClockEta(secondsFromNow) {
    const d = new Date(Date.now() + secondsFromNow * 1000);
    return d.toLocaleTimeString([], { hour: "2-digit", minute: "2-digit", second: "2-digit" });
  },

  parseWorkflow(line) {
    const m = line.match(/\[\s*(\d+)%\]\s+(\d+)\s+of\s+(\d+)/i);
    if (!m) {
      return false;
    }
    if (!this.workflow.startedAt) {
      this.workflow.startedAt = Date.now();
    }
    this.workflow.percent = Number(m[1]);
    this.workflow.completed = Number(m[2]);
    this.workflow.total = Number(m[3]);
    return true;
  },

  workflowEtaSeconds() {
    const wf = this.workflow;
    if (!wf.startedAt || wf.total <= 0) {
      return null;
    }
    const elapsed = (Date.now() - wf.startedAt) / 1000;
    if (wf.percent > 0 && wf.percent < 100) {
      return elapsed * (100 - wf.percent) / wf.percent;
    }
    if (wf.completed > 0 && wf.completed < wf.total) {
      return (elapsed / wf.completed) * (wf.total - wf.completed);
    }
    return null;
  },

  computeEta(proc) {
    if (proc.done || proc.percent >= 100 || (proc.total > 0 && proc.completed >= proc.total)) {
      return { label: "done", seconds: 0 };
    }

    const elapsed = (Date.now() - proc.startedAt) / 1000;

    if (proc.percent > 0 && proc.percent < 100) {
      return { label: "eta", seconds: elapsed * (100 - proc.percent) / proc.percent };
    }
    if (proc.completed > 0 && proc.total > proc.completed) {
      return {
        label: "eta",
        seconds: (elapsed / proc.completed) * (proc.total - proc.completed)
      };
    }

    const hint = this.hintFor(proc.name);
    if (hint) {
      const remaining = Math.max(8, hint - elapsed);
      return { label: "hint", seconds: remaining, elapsed };
    }

    const wfEta = this.workflowEtaSeconds();
    if (wfEta != null) {
      const active = [...this.processes.values()].filter((p) => !p.done).length || 1;
      return { label: "wf", seconds: wfEta / active, elapsed };
    }

    if (elapsed < 5) {
      return { label: "hint", seconds: hint, elapsed };
    }

    const remaining = Math.max(5, hint - elapsed);
    return { label: "hint", seconds: remaining, elapsed };
  },

  etaPrefix(proc) {
    if (proc.done) {
      const took = proc.frozenElapsed != null
        ? proc.frozenElapsed
        : (Date.now() - proc.startedAt) / 1000;
      return `done in ${this.formatDuration(took)}`;
    }

    const { label, seconds, elapsed } = this.computeEta(proc);
    const runFor = this.formatDuration(
      elapsed != null ? elapsed : (Date.now() - proc.startedAt) / 1000
    );

    if (seconds != null && seconds > 0) {
      const left = this.formatDuration(seconds);
      const clock = this.formatClockEta(seconds);
      return `~${left} left @ ${clock} · ran ${runFor}`;
    }

    return `~… left · ran ${runFor}`;
  },

  parse(line) {
    const t = line.trim();
    if (!t) {
      return null;
    }

    this.parseWorkflow(t);

    let m = t.match(
      /^\[([0-9a-f/]+)\]\s+process\s+>\s+(\S+)(?:\s+\(([^)]*)\))?\s+\[\s*(\d+)%\]\s+(\d+)\s+of\s+(\d+)(?:\s*✔)?/i
    );
    if (m) {
      const percent = Number(m[4]);
      const completed = Number(m[5]);
      const total = Number(m[6]);
      this.hashToName.set(m[1], m[2]);
      return {
        kind: "progress",
        id: m[1],
        name: m[2],
        tag: m[3] || "",
        percent,
        completed,
        total,
        done: percent >= 100 || (total > 0 && completed >= total) || /✔/.test(t),
        cached: /cached/i.test(t)
      };
    }

    m = t.match(/^\[([0-9a-f/]+)\].*\[\s*(\d+)%\]\s+(\d+)\s+of\s+(\d+)/i);
    if (m) {
      const name = this.hashToName.get(m[1]) || m[1];
      const percent = Number(m[2]);
      const completed = Number(m[3]);
      const total = Number(m[4]);
      return {
        kind: "progress",
        id: m[1],
        name,
        tag: "",
        percent,
        completed,
        total,
        done: percent >= 100 || (total > 0 && completed >= total) || /✔/.test(t),
        cached: /cached/i.test(t)
      };
    }

    m = t.match(/^\[([0-9a-f/]+)\]\s+(Submitted|Cached)\s+process\s+>\s+(\S+)\s*(?:\(([^)]*)\))?/i);
    if (m) {
      this.hashToName.set(m[1], m[3]);
      return {
        kind: m[2].toLowerCase(),
        id: m[1],
        name: m[3],
        tag: m[4] || "",
        percent: m[2].toLowerCase() === "cached" ? 100 : 0,
        completed: m[2].toLowerCase() === "cached" ? 1 : 0,
        total: 1,
        done: m[2].toLowerCase() === "cached",
        cached: m[2].toLowerCase() === "cached"
      };
    }

    m = t.match(/^\[([0-9a-f/]+)\]\s+process\s+>\s+(\S+)(?:\s+\(([^)]*)\))?/i);
    if (m) {
      this.hashToName.set(m[1], m[2]);
      return {
        kind: "started",
        id: m[1],
        name: m[2],
        tag: m[3] || "",
        percent: 0,
        completed: 0,
        total: 0,
        done: false
      };
    }

    m = t.match(/^\[-+\s*\]\s+process\s+>\s+(\S+)/i);
    if (m) {
      return {
        kind: "queued",
        id: null,
        name: m[1],
        tag: "",
        percent: 0,
        completed: 0,
        total: 0,
        done: false
      };
    }

    m = t.match(/^executor\s*>\s*(.+)$/i);
    if (m) {
      return {
        kind: "executor",
        id: "_executor",
        name: "executor",
        tag: m[1].trim(),
        percent: 0,
        completed: 0,
        total: 0,
        done: false
      };
    }

    return null;
  },

  procKey(name) {
    return this.normalizeName(name);
  },

  upsert(parsed) {
    const nameKey = this.procKey(parsed.name);
    let proc = this.processes.get(nameKey);
    if (!proc) {
      proc = {
        id: nameKey,
        name: parsed.name,
        tag: "",
        cached: false,
        startedAt: Date.now(),
        percent: 0,
        completed: 0,
        total: 0,
        done: false,
        row: null,
        lastLine: ""
      };
      this.processes.set(nameKey, proc);
    }

    if (parsed.id && parsed.id !== "_executor") {
      this.hashToProc.set(parsed.id, proc);
      this.hashToName.set(parsed.id, parsed.name);
    }

    if (parsed.tag) {
      proc.tag = parsed.tag;
    }
    if (parsed.kind === "cached" || parsed.cached) {
      proc.cached = true;
    }
    if (parsed.percent > 0) {
      proc.percent = Math.max(proc.percent, parsed.percent);
    }
    if (parsed.total > 0) {
      proc.total = parsed.total;
    }
    if (parsed.completed > 0) {
      proc.completed = parsed.completed;
    }
    if (parsed.done) {
      this.markDone(proc);
    }
    return proc;
  },

  formatDisplay(proc) {
    const tag = proc.tag ? ` (${proc.tag})` : "";
    if (proc.id === "_executor") {
      return `executor > ${proc.tag || "local"}`;
    }
    if (proc.done) {
      const cacheNote = proc.cached ? ", cached" : "";
      return `process > ${proc.name}${tag} [100%] 1 of 1${cacheNote} ✔`;
    }
    if (proc.percent > 0 || proc.total > 0) {
      const pct = String(proc.percent).padStart(3, " ");
      return `process > ${proc.name}${tag} [${pct}%] ${proc.completed} of ${proc.total}`;
    }
    return `process > ${proc.name}${tag} · waiting`;
  },

  renderRow(proc, className) {
    if (!proc.row) {
      return;
    }
    const eta = document.createElement("span");
    eta.className = "log-eta" + (proc.done ? " is-done" : "");
    if (proc.id === "_executor") {
      eta.textContent = "";
    } else {
      eta.textContent = this.etaPrefix(proc);
    }

    const text = document.createElement("span");
    text.className = "log-text";
    text.textContent = this.formatDisplay(proc);

    proc.row.replaceChildren(eta, document.createTextNode(eta.textContent ? " " : ""), text);
    proc.row.className = "log-line " + className;
  },

  refreshAll(terminal) {
    this.processes.forEach((proc) => {
      if (proc.row && !proc.done) {
        this.renderRow(proc, terminal.classify(proc.lastLine));
      }
    });
    if (this.allDone()) {
      this.stopTick();
    }
  },

  isProcessLine(line) {
    const t = String(line || "").trim();
    return /^\[-+\s*\]\s+process\s+>/i.test(t) ||
      /^\[[0-9a-f/]+\]\s+(?:Submitted|Cached\s+)?process\s+>/i.test(t) ||
      /^executor\s*>/i.test(t);
  },

  handleLine(terminal, line, className) {
    this.scanLineForCompletion(line, terminal);

    const parsed = this.parse(line);
    if (!parsed) {
      if (this.parseWorkflow(line)) {
        this.ensureTick(terminal);
        this.refreshAll(terminal);
      }
      return false;
    }

    this.ensureTick(terminal);
    const proc = this.upsert(parsed);
    proc.lastLine = line;

    if (!proc.row) {
      proc.row = document.createElement("div");
      proc.row.dataset.bfProcess = proc.id;
      terminal.logEl.appendChild(proc.row);
    }

    this.renderRow(proc, className);

    while (terminal.logEl.childNodes.length > terminal.maxLines) {
      terminal.logEl.removeChild(terminal.logEl.firstChild);
    }
    terminal.logEl.scrollTop = terminal.logEl.scrollHeight;

    if (proc.id === "_executor") {
      terminal.setStatus(`executor > ${proc.tag}`, "run");
    } else if (proc.done) {
      terminal.setStatus(`[${proc.name}] done in ${this.formatDuration(proc.frozenElapsed || 0)}`, "ok");
      if (this.allDone()) {
        this.stopTick();
      }
    } else {
      terminal.setStatus(`[${proc.name}] ${this.etaPrefix(proc)}`, "run");
    }

    return true;
  }
};

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

  isNextflowDebugLine(line) {
    return /\[Task monitor\]|TaskPollingMonitor|TaskHandler\[/i.test(line);
  },

  isRealError(line) {
    if (/error:\s*-/i.test(line) || this.isNextflowDebugLine(line)) {
      return false;
    }
    return /\b(ERROR|Error executing|terminated with an error|ModuleNotFoundError|Traceback|Command exit status|Command error)\b/i.test(line);
  },

  classify(line) {
    if (this.isNextflowDebugLine(line)) {
      return "log-muted";
    }
    if (this.isRealError(line)) {
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
    if (/Submitted process|Cached process|process >/.test(line) || (/\]\s+\S+/.test(line) && /\[\w{2}\/\w+\]/.test(line))) {
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
    if (this.isNextflowDebugLine(t)) {
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
    return /\b(ERROR|WARN)\b|✔|Launching|Process completed|BactFlow:|started|failed|ready in|Using Java|Using Nextflow|N E X T F L O W|executor >|Submitted process|Cached process|\[100%\]|\[[ ]*\d+%\]|\[  0%\]|Caused by:|Command exit status|running SPAdes|running Unicycler|running Flye|Running circlator|Stream disconnected|process >/i.test(t);
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
    BactflowProcessEta.reset();
    this.setStatus("Ready", "");
  },

  append(raw, force) {
    if (!this.logEl) {
      this.init();
    }
    const line = this.stripAnsi(raw);
    const trimmed = line.trim();

    if (trimmed === "Process completed" || /Nextflow finished successfully/i.test(trimmed)) {
      BactflowProcessEta.finalizeAll(this, "completed");
      if (trimmed !== "Process completed") {
        const row = document.createElement("div");
        row.className = "log-line log-ok";
        row.textContent = line;
        this.logEl.appendChild(row);
      }
      return;
    }

    BactflowProcessEta.scanLineForCompletion(line, this);

    const className = this.classify(line);
    if (BactflowProcessEta.handleLine(this, line, className)) {
      return;
    }

    if (!force && !this.isInteresting(line)) {
      BactflowProcessEta.parseWorkflow(line);
      return;
    }

    if (!force && BactflowProcessEta.isProcessLine(line)) {
      return;
    }

    const row = document.createElement("div");
    row.className = "log-line " + className;
    row.textContent = line || " ";
    this.logEl.appendChild(row);

    while (this.logEl.childNodes.length > this.maxLines) {
      this.logEl.removeChild(this.logEl.firstChild);
    }
    this.logEl.scrollTop = this.logEl.scrollHeight;

    if (/BactFlow: cleaning temporary/i.test(line)) {
      this.setStatus(line.replace(/\s+/g, " ").slice(0, 140), "ok");
    } else if (this.isRealError(line)) {
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

document.addEventListener("DOMContentLoaded", () => BactflowTerminal.init());
