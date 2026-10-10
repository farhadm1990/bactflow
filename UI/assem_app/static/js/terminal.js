const BactflowProcessEta = {
  processes: new Map(),
  hashToProc: new Map(),
  hashToName: new Map(),
  durations: new Map(),
  workflow: { startedAt: null, percent: 0, completed: 0, total: 0 },
  tickTimer: null,

  prettyName(name) {
    const labels = {
      envsetup: "Environment setup",
      testify: "Environment check",
      coveragefilt: "Coverage filter",
      nanoreadfilt: "Read filter",
      fastqconcater: "Read concat",
      deduper: "Deduper",
      assemblyflye1: "Flye",
      assemblyflye2: "Flye",
      assemblyspades: "SPAdes",
      assemblyunicycler: "Unicycler",
      assemblypacbio: "Flye (PacBio)",
      pacbioreadcheck: "PacBio read check",
      circulator: "Circulator",
      quastcheck: "QUAST",
      baktaannot: "Bakta",
      taxonomygtdbtk: "GTDB-Tk",
      checkmlineage: "CheckM",
      plasmiddetect: "geNomad"
    };
    const key = this.normalizeName(name);
    if (labels[key]) {
      return labels[key];
    }
    if (this.isHashId(name)) {
      return name;
    }
    return String(name || "process");
  },

  isHashId(value) {
    return /^[0-9a-f]{2,}\/[0-9a-f]+$/i.test(String(value || ""));
  },

  parseNameAndTag(raw) {
    const chunk = String(raw || "").trim();
    const tagged = chunk.match(/^([A-Za-z_]\w*)(?:\s+\((.*)\))?$/);
    if (tagged) {
      return { name: tagged[1], tag: tagged[2] || "" };
    }
    const ident = chunk.match(/([A-Za-z_]\w*)/);
    return { name: ident ? ident[1] : chunk, tag: "" };
  },

  hintsSec: {
    envsetup: 240,
    testify: 45,
    coverage_filt: 180,
    readfilt: 180,
    concat_reads: 120,
    assembly_flye1: 1500,
    assembly_flye2: 1500,
    assembly_pacbio: 1500,
    pacbio_read_check: 60,
    assembly_spades: 2100,
    assembly_unicycler: 2700,
    circulator: 300,
    quast: 240,
    quastcheck: 240,
    bakta: 900,
    gtdbtk: 600,
    checkm: 600,
    plasmiddetect: 900,
    genomad: 900
  },

  reset() {
    this.processes.clear();
    this.hashToProc.clear();
    this.hashToName.clear();
    this.durations.clear();
    this.loadLearnedDurations();
    this.workflow = { startedAt: null, percent: 0, completed: 0, total: 0 };
    if (this.tickTimer) {
      clearInterval(this.tickTimer);
      this.tickTimer = null;
    }
    this.clearPanel(typeof BactflowTerminal !== "undefined" ? BactflowTerminal : null);
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

  loadLearnedDurations() {
    try {
      const raw = localStorage.getItem("bactflow_process_durations_v1");
      if (!raw) {
        return;
      }
      const obj = JSON.parse(raw);
      Object.entries(obj || {}).forEach(([k, v]) => {
        const sec = Number(v);
        if (k && Number.isFinite(sec) && sec >= 20) {
          this.durations.set(k, sec);
        }
      });
    } catch (e) { /* ignore */ }
  },

  saveLearnedDurations() {
    try {
      const obj = {};
      this.durations.forEach((v, k) => {
        obj[k] = Math.round(v);
      });
      localStorage.setItem("bactflow_process_durations_v1", JSON.stringify(obj));
    } catch (e) { /* ignore */ }
  },

  ensurePanel(terminal) {
    if (!terminal) {
      return;
    }
    if (!terminal.processPanelEl) {
      terminal.processPanelEl = document.getElementById("bf-process-panel");
    }
    if (!terminal.processListEl) {
      terminal.processListEl = document.getElementById("bf-process-list");
    }
    if (!terminal.processOverallFillEl) {
      terminal.processOverallFillEl = document.getElementById("bf-process-overall-fill");
    }
    if (!terminal.processOverallMetaEl) {
      terminal.processOverallMetaEl = document.getElementById("bf-process-overall-meta");
    }
    if (terminal.processPanelEl) {
      terminal.processPanelEl.style.display = "block";
    }
  },

  clearPanel(terminal) {
    const panel = (terminal && terminal.processPanelEl) || document.getElementById("bf-process-panel");
    const list = (terminal && terminal.processListEl) || document.getElementById("bf-process-list");
    const fill = (terminal && terminal.processOverallFillEl) || document.getElementById("bf-process-overall-fill");
    const meta = (terminal && terminal.processOverallMetaEl) || document.getElementById("bf-process-overall-meta");
    if (list) {
      list.replaceChildren();
    }
    if (fill) {
      fill.style.width = "0%";
      fill.className = "bf-process-fill is-workflow";
    }
    if (meta) {
      meta.textContent = "—";
    }
    if (panel) {
      panel.style.display = "none";
    }
  },

  displayPercent(proc) {
    if (!proc || proc.id === "_executor") {
      return 0;
    }
    if (proc.done || proc.cached) {
      return 100;
    }
    if (this.isWaiting(proc)) {
      return 0;
    }
    if (proc.percent > 0 && proc.percent < 100) {
      return Math.max(1, Math.min(99, proc.percent));
    }
    const elapsed = this.elapsedSeconds(proc);
    const expected = Math.max(30, this.hintFor(proc.name) || 180);
    let p = 95 * (1 - Math.exp(-elapsed / expected));
    if (proc.total > 0) {
      // "0 of 1" stays put while a single task runs, so never let the task count drag the estimate down.
      p = Math.max(p, (Math.max(0, proc.completed) / proc.total) * 100);
    }
    return Math.max(1, Math.min(99, Math.round(p)));
  },

  remainingSeconds(proc) {
    if (!proc || proc.done) {
      return 0;
    }
    if (this.isWaiting(proc)) {
      return null;
    }
    const pct = this.displayPercent(proc);
    const elapsed = this.elapsedSeconds(proc);
    if (pct >= 5 && pct < 100) {
      return Math.max(1, elapsed * (100 - pct) / pct);
    }
    const expected = Math.max(30, this.hintFor(proc.name) || 180);
    return Math.max(1, expected - elapsed);
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
    if (!proc.cached && proc.frozenElapsed >= 45) {
      this.durations.set(this.normalizeName(proc.name), proc.frozenElapsed);
      this.saveLearnedDurations();
    }
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

    const elapsed = this.elapsedSeconds(proc);

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
    if (hint && elapsed < hint) {
      return { label: "hint", seconds: hint - elapsed, elapsed };
    }
    if (hint && elapsed >= hint) {
      const extra = Math.max(20, elapsed * 0.2);
      return { label: "overrun", seconds: extra, elapsed };
    }

    const wfEta = this.workflowEtaSeconds();
    if (wfEta != null && wfEta > 0) {
      const active = [...this.processes.values()].filter((p) => !p.done && p.id !== "_executor").length || 1;
      return { label: "wf", seconds: wfEta / active, elapsed };
    }

    return { label: "running", seconds: null, elapsed };
  },

  elapsedSeconds(proc) {
    if (proc.done && proc.frozenElapsed != null) {
      return proc.frozenElapsed;
    }
    const start = proc.runningAt || proc.startedAt;
    return (Date.now() - start) / 1000;
  },

  isWaiting(proc) {
    return !proc.done && proc.phase !== "running";
  },

  etaPrefix(proc) {
    const runFor = this.formatDuration(this.elapsedSeconds(proc));
    if (proc.done) {
      return `done in ${runFor}`;
    }
    if (this.isWaiting(proc)) {
      return `waiting · ${runFor}`;
    }
    const rem = this.remainingSeconds(proc);
    const pct = this.displayPercent(proc);
    if (rem != null && rem > 0) {
      return `${pct}% · ${runFor} · ETA ${this.formatDuration(rem)} (~${this.formatClockEta(rem)})`;
    }
    return `${pct}% · ${runFor}`;
  },

  parse(line) {
    const t = line.trim();
    if (!t) {
      return null;
    }

    this.parseWorkflow(t);

    const progress = (hash, name, tag, percent, completed, total) => {
      const resolved = (!this.isHashId(name) && name)
        ? name
        : (this.hashToName.get(hash) || name);
      if (hash && resolved && !this.isHashId(resolved)) {
        this.hashToName.set(hash, resolved);
      }
      return {
        kind: "progress",
        id: hash,
        name: resolved,
        tag: tag || "",
        percent,
        completed,
        total,
        done: percent >= 100 || (total > 0 && completed >= total) || /✔/.test(t),
        cached: /cached/i.test(t)
      };
    };

    let m = t.match(
      /^\[([0-9a-f/]+)\]\s+process\s+>\s+(\S+)(?:\s+\(([^)]*)\))?\s+\[\s*(\d+)%\]\s+(\d+)\s+of\s+(\d+)/i
    );
    if (m) {
      return progress(m[1], m[2], m[3], Number(m[4]), Number(m[5]), Number(m[6]));
    }

    m = t.match(
      /^\[([0-9a-f/]+)\]\s+(.+?)\s+\[\s*(\d+)%\]\s+(\d+)\s+of\s+(\d+)/i
    );
    if (m) {
      const middle = m[2].trim().replace(/^process\s*>\s*/i, "");
      const { name, tag } = this.parseNameAndTag(middle);
      return progress(m[1], name, tag, Number(m[3]), Number(m[4]), Number(m[5]));
    }

    // Nextflow's live (ANSI) view: "[ec/570f4f] testify | 1 of 1 ✔"; long names are shortened with "…".
    m = t.match(/^\[([0-9a-f]{2}\/[0-9a-f]{6})\]\s+(.+?)\s+\|\s+(\d+)\s+of\s+(\d+)/i);
    if (m) {
      let name = this.hashToName.get(m[1]);
      if (!name && !m[2].includes("\u2026")) {
        name = this.parseNameAndTag(m[2].replace(/^process\s*>\s*/i, "")).name;
      }
      if (!name) {
        const prefix = m[2].split("\u2026")[0].trim().toLowerCase();
        const hits = [...this.processes.values()].filter((p) =>
          prefix && p.id !== "_executor" && (!p.workHash || p.workHash === m[1]) &&
          String(p.name).toLowerCase().startsWith(prefix)
        );
        if (hits.length === 1) {
          name = hits[0].name;
        }
      }
      if (!name) {
        return null;
      }
      const completed = Number(m[3]);
      const total = Number(m[4]);
      const done = /✔/.test(t) || (total > 0 && completed >= total);
      return progress(m[1], name, "", done ? 100 : Math.floor(100 * completed / Math.max(1, total)), completed, total);
    }

    m = t.match(/\[([0-9a-f/]+)\]\s+(Submitted|Cached)\s+process\s+>\s+(\S+)\s*(?:\(([^)]*)\))?/i);
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

    m = t.match(/\[([0-9a-f/]+)\]\s+process\s+>\s+(\S+)(?:\s+\(([^)]*)\))?/i);
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

    m = t.match(/^\[-+\s*\]\s+process\s+>\s+(\S+)(?:\s+\(([^)]*)\))?/i);
    if (m) {
      return {
        kind: "queued",
        id: null,
        name: m[1],
        tag: m[2] || "",
        percent: 0,
        completed: 0,
        total: 0,
        done: false
      };
    }

    m = t.match(/^\[-+\s*\]\s+([A-Za-z_]\w*)(?:\s+\(([^)]*)\))?(?:\s+-+\s*)?$/);
    if (m && !/^executor$/i.test(m[1])) {
      return {
        kind: "queued",
        id: null,
        name: m[1],
        tag: m[2] || "",
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

    m = t.match(/^running\s+(circlator|spades|unicycler|flye)\b/i);
    if (m) {
      const tool = m[1].toLowerCase();
      let name = tool;
      if (tool === "spades") {
        name = "assembly_spades";
      } else if (tool === "unicycler") {
        name = "assembly_unicycler";
      } else if (tool === "flye") {
        name = /pacbio/i.test(t) ? "assembly_pacbio" : "assembly_flye1";
      } else {
        name = "circulator";
      }
      return {
        kind: "started",
        id: null,
        name,
        tag: "",
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

  displayName(proc) {
    if (!proc || proc.id === "_executor") {
      return "executor";
    }
    const pretty = this.prettyName(proc.name);
    const hash = proc.workHash || "";
    return hash ? `${pretty} [${hash}]` : pretty;
  },

  upsert(parsed) {
    if (this.isHashId(parsed.name) && parsed.id && this.hashToName.has(parsed.id)) {
      parsed.name = this.hashToName.get(parsed.id);
    }

    let proc = (parsed.id && parsed.id !== "_executor")
      ? this.hashToProc.get(parsed.id)
      : null;
    const nameKey = this.procKey(parsed.name);
    if (!proc && !this.isHashId(parsed.name)) {
      proc = this.processes.get(nameKey);
    }
    if (!proc) {
      proc = {
        id: nameKey,
        name: parsed.name,
        tag: "",
        cached: false,
        startedAt: Date.now(),
        runningAt: null,
        phase: "waiting",
        percent: 0,
        completed: 0,
        total: 0,
        done: false,
        row: null,
        lastLine: ""
      };
      this.processes.set(nameKey, proc);
    }

    if (parsed.name && !this.isHashId(parsed.name)) {
      proc.name = parsed.name;
      const betterKey = this.procKey(parsed.name);
      if (proc.id !== betterKey) {
        this.processes.delete(proc.id);
        proc.id = betterKey;
        this.processes.set(betterKey, proc);
      }
    }

    if (parsed.id && parsed.id !== "_executor") {
      this.hashToProc.set(parsed.id, proc);
      if (this.isHashId(parsed.id)) {
        proc.workHash = parsed.id;
      }
      if (!this.isHashId(parsed.name)) {
        this.hashToName.set(parsed.id, parsed.name);
      }
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
    if (parsed.kind === "started" || parsed.kind === "progress" || parsed.kind === "submitted") {
      proc.phase = "running";
      if (!proc.runningAt) {
        proc.runningAt = Date.now();
      }
    } else if (parsed.kind === "queued") {
      if (proc.phase !== "running") {
        proc.phase = "waiting";
      }
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
    const title = this.displayName(proc) + tag;
    if (proc.done) {
      const cacheNote = proc.cached ? " · cached" : "";
      return `${title}${cacheNote} ✔`;
    }
    if (this.isWaiting(proc)) {
      return `${title} · waiting`;
    }
    if (proc.percent > 0 || proc.total > 0) {
      const pct = String(proc.percent).padStart(3, " ");
      return `${title} [${pct}%] ${proc.completed} of ${proc.total}`;
    }
    return `${title} · running`;
  },

  renderRow(proc, className) {
    if (!proc || proc.id === "_executor") {
      return;
    }
    const terminal = typeof BactflowTerminal !== "undefined" ? BactflowTerminal : null;
    this.ensurePanel(terminal);
    const list = terminal && terminal.processListEl;
    if (!list) {
      return;
    }
    if (!proc.row || !proc.row.isConnected) {
      proc.row = document.createElement("div");
      proc.row.className = "bf-proc-card";
      proc.row.dataset.bfProcess = proc.id;
      list.appendChild(proc.row);
    }

    const pct = this.displayPercent(proc);
    const waiting = this.isWaiting(proc);
    const state = proc.done ? (proc.cached ? "cached" : "done") : waiting ? "wait" : "run";
    proc.row.className = "bf-proc-card" + (state === "done" || state === "cached" ? " is-done" : state === "run" ? " is-run" : "");

    const tag = proc.tag ? ` (${proc.tag})` : "";
    const title = this.prettyName(proc.name) + tag;
    const badgeLabel = state === "cached" ? "cached" : state === "done" ? "done" : state === "wait" ? "queued" : "running";
    let meta = "";
    if (proc.done) {
      meta = `finished in ${this.formatDuration(proc.frozenElapsed || 0)}`;
    } else if (waiting) {
      meta = `queued · ${this.formatDuration(this.elapsedSeconds(proc))}`;
    } else {
      const rem = this.remainingSeconds(proc);
      const runFor = this.formatDuration(this.elapsedSeconds(proc));
      meta = rem != null
        ? `elapsed ${runFor} · ETA ${this.formatDuration(rem)} (~${this.formatClockEta(rem)})`
        : `elapsed ${runFor}`;
      if (proc.total > 0) {
        meta = `${proc.completed}/${proc.total} · ` + meta;
      }
    }

    proc.row.innerHTML = `
      <div class="bf-proc-head">
        <div>
          <span class="bf-proc-name"></span>
          <span class="bf-proc-badge is-${state}"></span>
        </div>
        <span class="bf-proc-pct"></span>
      </div>
      <div class="bf-process-track"><div class="bf-process-fill is-${state}"></div></div>
      <div class="bf-proc-meta" style="margin-top:5px;"></div>
    `;
    proc.row.querySelector(".bf-proc-name").textContent = title;
    proc.row.querySelector(".bf-proc-badge").textContent = badgeLabel;
    proc.row.querySelector(".bf-proc-pct").textContent = `${pct}%`;
    proc.row.querySelector(".bf-process-fill").style.width = (waiting ? 100 : pct) + "%";
    proc.row.querySelector(".bf-proc-meta").textContent = meta;
    this.renderWorkflowBar(terminal);
  },

  renderWorkflowBar(terminal) {
    if (!terminal) {
      return;
    }
    this.ensurePanel(terminal);
    const fill = terminal.processOverallFillEl;
    const meta = terminal.processOverallMetaEl;
    if (!fill || !meta) {
      return;
    }
    const procs = [...this.processes.values()].filter((p) => p.id !== "_executor");
    let pct = 0;
    if (this.workflow.total > 0) {
      pct = Math.max(0, Math.min(100, this.workflow.percent || Math.round(100 * this.workflow.completed / this.workflow.total)));
      meta.textContent = `${this.workflow.completed}/${this.workflow.total} · ${pct}%`;
    } else if (procs.length) {
      const sum = procs.reduce((acc, p) => acc + this.displayPercent(p), 0);
      pct = Math.round(sum / procs.length);
      const done = procs.filter((p) => p.done).length;
      meta.textContent = `${done}/${procs.length} tasks · ${pct}%`;
    } else {
      meta.textContent = "—";
    }
    fill.style.width = pct + "%";
    fill.className = "bf-process-fill is-workflow" + (pct >= 100 ? " is-done" : "");
  },

  refreshAll(terminal) {
    this.processes.forEach((proc) => {
      if (proc.row && proc.id !== "_executor") {
        this.renderRow(proc, terminal ? terminal.classify(proc.lastLine) : "log-process");
      }
    });
    this.renderWorkflowBar(terminal);
    if (this.allDone()) {
      this.stopTick();
    }
  },

  applyLiveActivity(line, terminal) {
    const t = String(line || "");
    const rules = [
      [/running SPAdes|spades.py/i, /spades/i],
      [/running Unicycler|unicycler/i, /unicycler/i],
      [/running Flye|flye /i, /flye|pacbio/i],
      [/Running circlator|circlator fixstart/i, /circul/i],
      [/quast\.py|Running QUAST/i, /quast/i],
      [/geNomad|plasmidDetect|plasmid_detect/i, /plasmid|genomad/i],
      [/checkm /i, /checkm/i],
      [/gtdbtk|GTDB-Tk/i, /gtdb|taxonomy/i],
      [/bakta /i, /bakta/i]
    ];
    let hit = false;
    for (const [re, nameRe] of rules) {
      if (!re.test(t)) {
        continue;
      }
      this.processes.forEach((proc) => {
        if (proc.done || proc.id === "_executor") {
          return;
        }
        if (nameRe.test(`${proc.name} ${proc.id}`)) {
          proc.phase = "running";
          if (!proc.runningAt) {
            proc.runningAt = Date.now();
          }
          proc.lastLine = line;
          if (proc.row && terminal) {
            this.renderRow(proc, terminal.classify(line));
          }
          hit = true;
        }
      });
    }
    return hit;
  },

  isProcessLine(line) {
    const t = String(line || "").trim();
    return /^\[-+\s*\]\s+(?:process\s+>\s+)?[A-Za-z_]/i.test(t) ||
      /^\[[0-9a-f/]+\]\s+(?:Submitted|Cached\s+)?(?:process\s+>\s+)?[A-Za-z_]/i.test(t) ||
      /^executor\s*>/i.test(t) ||
      /^running\s+(circlator|spades|unicycler|flye)\b/i.test(t);
  },

  handleLine(terminal, line, className) {
    this.scanLineForCompletion(line, terminal);
    this.applyLiveActivity(line, terminal);

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

    if (proc.id !== "_executor") {
      this.renderRow(proc, className);
    } else {
      this.renderWorkflowBar(terminal);
    }

    const shown = this.prettyName(proc.name);
    if (proc.id === "_executor") {
      terminal.setStatus(`executor > ${proc.tag}`, "run");
    } else if (proc.done) {
      terminal.setStatus(`${shown} · done in ${this.formatDuration(proc.frozenElapsed || 0)}`, "ok");
      if (this.allDone()) {
        this.stopTick();
        this.renderWorkflowBar(terminal);
      }
    } else if (this.isWaiting(proc)) {
      terminal.setStatus(`${shown} · queued`, "run");
    } else {
      terminal.setStatus(`${shown} · ${this.displayPercent(proc)}%`, "run");
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
    this.processPanelEl = document.getElementById("bf-process-panel");
    this.processListEl = document.getElementById("bf-process-list");
    this.processOverallFillEl = document.getElementById("bf-process-overall-fill");
    this.processOverallMetaEl = document.getElementById("bf-process-overall-meta");
    BactflowProcessEta.loadLearnedDurations();
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
    if (/^\[-\s*\]/.test(t) && !/process\s*>/i.test(t)) {
      return false;
    }
    if (/Tip: you can|Check '\.nextflow\.log'|WORKFLOW OUTPUT DEFINITION|is available - Please consider/.test(t)) {
      return false;
    }
    return /\b(ERROR|WARN)\b|✔|Launching|Process completed|BactFlow:|started|failed|ready in|Using Java|Using Nextflow|N E X T F L O W|executor >|Submitted process|Cached process|\[100%\]|\[[ ]*\d+%\]|\[  0%\]|Caused by:|Command exit status|running SPAdes|running Unicycler|running Flye|Running circlator|geNomad|plasmidDetect|Stream disconnected|process >|\[[0-9a-f]{2}\/[0-9a-f]+\]/i.test(t);
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
    BactflowProcessEta.clearPanel(this);
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

  attachStream(eventSource) {
    if (!eventSource || eventSource._bfMetersBound) {
      return;
    }
    eventSource._bfMetersBound = true;
    eventSource.addEventListener("stats", (ev) => {
      try {
        this.apply(JSON.parse(ev.data));
      } catch (err) {
        // ignore malformed stats frames
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

  apply(s) {
    if (!s || s.ram_percent == null && s.cpu_percent == null) {
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

  markWaiting(el, msg) {
    if (!el) {
      return;
    }
    el.classList.add("is-missing");
    el.textContent = msg || "waiting…";
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
      for (const url of [
        "/sys_stats",
        "/bactflow_status",
        `/static/sys_stats.json?t=${stamp}`
      ]) {
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
  }
};

window.BactflowMeters = BactflowMeters;

document.addEventListener("DOMContentLoaded", () => {
  BactflowTerminal.init();
  BactflowMeters.init();
});
