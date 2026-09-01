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
    if (/Submitted process|Cached process|Submitted process/.test(line) || /\]\s+\S+/.test(line) && /\[\w{2}\/\w+\]/.test(line)) {
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
    return /ERROR|WARN|✔|Launching|Process completed|started|failed|ready in|Using Java|N E X T F L O W|executor >|Submitted process|Cached process|\[100%\]|\[  0%\]|Caused by:|Command exit status|running SPAdes|running Unicycler|running Flye|Running circlator|Stream disconnected/i.test(t);
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

document.addEventListener("DOMContentLoaded", () => BactflowTerminal.init());
