/* BactFlow Run Dashboard — option-driven pipeline, live run state, resources, theme */
(function () {
  const $ = (sel, root) => (root || document).querySelector(sel);

  const ICONS = {
    setup:
      '<svg viewBox="0 0 24 24"><path d="M14.7 6.3a4 4 0 0 0-5.4 5.4L4 17v3h3l5.3-5.3a4 4 0 0 0 5.4-5.4l-2.5 2.5-2.4-.6-.6-2.4z"/></svg>',
    check:
      '<svg viewBox="0 0 24 24"><path d="M12 3l7 3v5c0 4.5-3 8.3-7 10-4-1.7-7-5.5-7-10V6z"/><path d="M9 12l2 2 4-4"/></svg>',
    pool:
      '<svg viewBox="0 0 24 24"><path d="M4 6h16M4 12h16M4 18h16"/><path d="M8 3v6M16 9v6M11 15v6"/></svg>',
    dedup:
      '<svg viewBox="0 0 24 24"><rect x="4" y="4" width="11" height="11" rx="2"/><path d="M9 20h9a2 2 0 0 0 2-2V9"/></svg>',
    filter:
      '<svg viewBox="0 0 24 24"><path d="M4 5h16l-6 7v6l-4 2v-8z"/></svg>',
    coverage:
      '<svg viewBox="0 0 24 24"><path d="M3 18c3-8 5-12 9-12s6 4 9 12"/><path d="M3 18h18"/></svg>',
    readcheck:
      '<svg viewBox="0 0 24 24"><circle cx="11" cy="11" r="6"/><path d="M20 20l-4.5-4.5M8.5 11l2 2 3-3"/></svg>',
    flye:
      '<svg viewBox="0 0 24 24"><path d="M3 12c3-5 6-5 9 0s6 5 9 0"/><path d="M3 17c3-5 6-5 9 0s6 5 9 0" opacity=".5"/><path d="M3 7c3-5 6-5 9 0s6 5 9 0" opacity=".5"/></svg>',
    spades:
      '<svg viewBox="0 0 24 24"><circle cx="6" cy="6" r="2.5"/><circle cx="18" cy="6" r="2.5"/><circle cx="12" cy="18" r="2.5"/><path d="M8 7.5l3 8M16 7.5l-3 8M8.5 6h7"/></svg>',
    unicycler:
      '<svg viewBox="0 0 24 24"><circle cx="12" cy="12" r="8"/><path d="M12 4v16M4 12h16"/></svg>',
    pacbio:
      '<svg viewBox="0 0 24 24"><path d="M4 12h16"/><path d="M7 8l-3 4 3 4M17 8l3 4-3 4"/></svg>',
    circulator:
      '<svg viewBox="0 0 24 24"><path d="M3 12a9 9 0 1 0 3-6.7"/><path d="M3 4v5h5"/></svg>',
    bakta:
      '<svg viewBox="0 0 24 24"><path d="M5 4h10l4 4v12H5z"/><path d="M8 11h8M8 15h6"/></svg>',
    gtdb:
      '<svg viewBox="0 0 24 24"><path d="M12 21v-7M12 14l-5-4M12 14l5-4M7 10V6M17 10V6M7 6l-2-2M7 6l2-2M17 6l-2-2M17 6l2-2"/></svg>',
    checkm:
      '<svg viewBox="0 0 24 24"><circle cx="12" cy="12" r="9"/><path d="M8 12l3 3 5-6"/></svg>',
    genomad:
      '<svg viewBox="0 0 24 24"><circle cx="12" cy="12" r="8"/><circle cx="12" cy="12" r="3"/><path d="M12 4v3M12 17v3"/></svg>',
    quast:
      '<svg viewBox="0 0 24 24"><path d="M4 19V5M4 19h16M8 15v4M12 11v8M16 8v11"/></svg>',
    done: '<svg viewBox="0 0 24 24"><path d="M5 12l5 5L20 7"/></svg>',
    fail: '<svg viewBox="0 0 24 24"><path d="M6 6l12 12M18 6L6 18"/></svg>',
  };

  const THEME_ICONS = {
    dark: '<svg viewBox="0 0 24 24"><path d="M20 14.5A8 8 0 0 1 9.5 4a8 8 0 1 0 10.5 10.5z"/></svg>',
    light:
      '<svg viewBox="0 0 24 24"><circle cx="12" cy="12" r="4"/><path d="M12 2v2M12 20v2M4.9 4.9l1.4 1.4M17.7 17.7l1.4 1.4M2 12h2M20 12h2M4.9 19.1l1.4-1.4M17.7 6.3l1.4-1.4"/></svg>',
  };

  const FAIL_RE =
    /ERROR ~|history\.lock|Execution aborted|terminated with an error|Error executing process|Failed to start|WorkflowScriptError|does not exist|No ONT FASTQ|No PacBio FASTQ|Nextflow exited with code|Process failed|Command exit status:\s*[1-9]|Command error|exited with code\s+[1-9]/i;

  function cleanFailMessage(line) {
    return String(line || "")
      .replace(/^ERROR\s*~\s*/i, "")
      .replace(/^BactFlow:\s*/i, "")
      .replace(/\s+/g, " ")
      .trim()
      .slice(0, 500);
  }

  const STEP_TARGETS = {
    envsetup: [".rd-logs"],
    testify: [".rd-logs"],
    concat: [".rd-logs"],
    pbcheck: [".rd-logs"],
    dedup: [".rd-logs"],
    readfilt: [".rd-logs"],
    covfilt: [".rd-logs"],
    asm: [".rd-logs"],
    circ: ["circFastaDiv", "circDiv", ".rd-logs"],
    bakta: ["baktaDiv", "#rd-results"],
    gtdb: ["taxa_class", "gtdbTreeDiv", "#rd-results"],
    checkm: ["checkmDiv", "checkmTreeDiv", "#rd-results"],
    genomad: ["plasmidDiv", "#rd-results"],
    quast: ["quastDiv", "#rd-results"],
  };

  function tracker() {
    return typeof BactflowProcessEta !== "undefined" ? BactflowProcessEta : null;
  }

  function terminal() {
    return typeof BactflowTerminal !== "undefined" ? BactflowTerminal : null;
  }

  function norm(s) {
    return String(s || "").toLowerCase().replace(/[^a-z0-9]+/g, "");
  }

  const RunDashboard = {
    startedAt: null,
    launchAt: null,
    endedAt: null,
    nfFinished: false,
    failed: false,
    failedProc: "",
    failMessage: "",
    doneAtFail: null,
    stopping: false,
    frozenRows: null,
    serverRunning: false,
    lastStats: null,
    plannedKey: null,
    autoScroll: true,

    init() {
      this.initTheme();
      this.bindUserMenu();
      this.bindLogTools();
      this.bindRunHooks();
      this.hookTerminal();
      this.hookMeters();
      this.bindPlotResize();
      this.bindStickyPipeline();
      this.bindLogsAlign();
      this.bindStepJumps();
      this.stampRunId();
      const form = this.form();
      if (form) {
        form.addEventListener("change", () => this.render());
        form.addEventListener("input", () => this.render());
      }
      setInterval(() => this.render(), 1000);
      this.render();
    },

    form() {
      return document.getElementById("runForm") || document.getElementById("postForm");
    },

    /* —— Theme —— */
    initTheme() {
      const btn = document.getElementById("rd-theme-bt");
      const apply = (theme) => {
        document.documentElement.setAttribute("data-theme", theme);
        if (btn) {
          const next = theme === "light" ? "dark" : "light";
          btn.innerHTML = THEME_ICONS[next] + (next === "light" ? "Bright" : "Dark");
          btn.title = "Switch to " + (next === "light" ? "bright" : "dark") + " theme";
        }
      };
      const order = ["dark", "plain", "light"];
      const labels = { dark: "Plasma", plain: "Plain dark", light: "Bright" };
      let theme = "dark";
      try {
        const saved = localStorage.getItem("bactflow_theme");
        theme = order.includes(saved) ? saved : "dark";
      } catch (e) { /* ignore */ }
      const paint = (value) => {
        document.documentElement.setAttribute("data-theme", value);
        if (btn) {
          const next = order[(order.indexOf(value) + 1) % order.length];
          btn.textContent = labels[next];
          btn.title = "Switch to " + labels[next];
        }
      };
      paint(theme);
      if (btn && btn.dataset.themeBound !== "1") {
        btn.dataset.themeBound = "1";
        btn.addEventListener("click", () => {
          theme = order[(order.indexOf(theme) + 1) % order.length];
          try {
            localStorage.setItem("bactflow_theme", theme);
          } catch (e) { /* ignore */ }
          paint(theme);
        });
      }
      this.bindFocus();
    },

    bindFocus() {
      const btn = document.getElementById("rd-focus-bt");
      const app = document.querySelector(".rd-app");
      if (!btn || !app || btn.dataset.focusBound === "1") {
        return;
      }
      btn.dataset.focusBound = "1";
      const apply = (on) => {
        app.classList.toggle("is-focus", on);
        btn.textContent = on ? "Show panels" : "Focus";
        btn.title = on ? "Show the setup menu and resource panels" : "Hide the setup menu and resource panels";
        if (this.alignLogsWithNextflow) {
          requestAnimationFrame(() => this.alignLogsWithNextflow());
        }
      };
      let on = false;
      try {
        on = localStorage.getItem("bactflow_focus") === "1";
      } catch (e) { /* ignore */ }
      apply(on);
      btn.addEventListener("click", () => {
        on = !app.classList.contains("is-focus");
        try {
          localStorage.setItem("bactflow_focus", on ? "1" : "0");
        } catch (e) { /* ignore */ }
        apply(on);
      });
    },

    /* —— Operator popup —— */
    bindUserMenu() {
      const btn = document.getElementById("rd-user-btn");
      const menu = document.getElementById("rd-user-menu");
      if (!btn || !menu || btn.dataset.menuBound === "1") {
        return;
      }
      btn.dataset.menuBound = "1";
      const place = () => {
        const r = btn.getBoundingClientRect();
        menu.style.position = "fixed";
        menu.style.left = Math.max(8, r.left) + "px";
        menu.style.width = Math.max(180, r.width) + "px";
        menu.style.right = "auto";
        menu.style.bottom = Math.max(8, window.innerHeight - r.top + 6) + "px";
        menu.style.top = "auto";
        menu.style.zIndex = "4000";
      };
      const setOpen = (open) => {
        menu.classList.toggle("is-open", open);
        btn.setAttribute("aria-expanded", open ? "true" : "false");
        if (open) {
          place();
        }
      };
      btn.addEventListener("click", (e) => {
        e.preventDefault();
        e.stopPropagation();
        setOpen(!menu.classList.contains("is-open"));
      });
      document.addEventListener("click", (e) => {
        if (!menu.contains(e.target) && !btn.contains(e.target)) {
          setOpen(false);
        }
      });
      document.addEventListener("keydown", (e) => {
        if (e.key === "Escape") {
          setOpen(false);
        }
      });
      window.addEventListener("resize", () => {
        if (menu.classList.contains("is-open")) {
          place();
        }
      });
    },

    /* Plotly sizes charts when they are drawn; charts drawn while hidden overflow until resized. */
    bindPlotResize() {
      const pane = document.getElementById("assem-results-pane");
      if (!pane) {
        return;
      }
      let timer = null;
      const resizeAll = () => {
        timer = null;
        if (!window.Plotly) {
          return;
        }
        pane.querySelectorAll(".js-plotly-plot").forEach((el) => {
          const w = el.clientWidth;
          const drawn = el._fullLayout && el._fullLayout.width;
          if (w && drawn && Math.abs(drawn - w) > 4) {
            try {
              Plotly.Plots.resize(el);
            } catch (e) { /* ignore */ }
          }
        });
      };
      this.schedulePlotResize = () => {
        if (!timer) {
          timer = setTimeout(resizeAll, 150);
        }
      };
      if (window.ResizeObserver) {
        new ResizeObserver(this.schedulePlotResize).observe(pane);
      }
      new MutationObserver(this.schedulePlotResize).observe(pane, {
        subtree: true,
        childList: true,
        attributes: true,
        attributeFilter: ["style", "class"],
      });
    },

    /* —— Log tools —— */
    bindLogTools() {
      const auto = $("#rd-autoscroll");
      if (auto) {
        auto.checked = true;
        auto.addEventListener("change", () => {
          this.autoScroll = !!auto.checked;
          if (window.BactflowTerminal) {
            window.BactflowTerminal.followScroll = this.autoScroll;
            if (this.autoScroll && typeof window.BactflowTerminal.scrollToBottom === "function") {
              window.BactflowTerminal.scrollToBottom();
            }
          }
        });
      }
      const dl = $("#rd-download-log");
      if (dl) {
        dl.addEventListener("click", () => {
          const body = document.getElementById("output-bactflow");
          const text = body ? body.innerText || body.textContent || "" : "";
          const a = document.createElement("a");
          a.href = URL.createObjectURL(new Blob([text], { type: "text/plain" }));
          a.download = "bactflow-log.txt";
          a.click();
          URL.revokeObjectURL(a.href);
        });
      }
    },

    scrollTo(sel) {
      const el = typeof sel === "string"
        ? (sel.charAt(0) === "#" || sel.charAt(0) === "." ? $(sel) : document.getElementById(sel))
        : sel;
      if (!el) {
        return;
      }
      const scroller = document.querySelector(".rd-center");
      const flow = document.getElementById("rd-flow-card");
      if (scroller && scroller.contains(el) && scroller.scrollHeight > scroller.clientHeight + 8) {
        if (flow && scroller.scrollTop > 8) {
          flow.classList.add("is-compact");
        }
        const offset = (flow ? flow.getBoundingClientRect().height : 0) + 10;
        const top =
          el.getBoundingClientRect().top -
          scroller.getBoundingClientRect().top +
          scroller.scrollTop -
          offset;
        scroller.scrollTo({ top: Math.max(0, top), behavior: "smooth" });
        return;
      }
      el.scrollIntoView({ behavior: "smooth", block: "start" });
    },

    bindStickyPipeline() {
      const scroller = document.querySelector(".rd-center");
      const card = document.getElementById("rd-flow-card");
      if (!scroller || !card) {
        return;
      }
      const update = () => {
        const compact = scroller.scrollTop > 16;
        if (card.classList.contains("is-compact") === compact) {
          return;
        }
        card.classList.toggle("is-compact", compact);
        requestAnimationFrame(() => this.alignLogsWithNextflow());
      };
      scroller.addEventListener("scroll", update, { passive: true });
      update();
    },

    bindStepJumps() {
      const onJump = (e) => {
        const hit = e.target.closest(".rd-stage-node, .rd-stage-name, .rd-progress-row");
        if (!hit) {
          return;
        }
        const el = hit.closest("[data-step]");
        if (!el || !el.classList.contains("is-jump")) {
          return;
        }
        e.preventDefault();
        this.jumpToStep(el.dataset.step);
      };
      const track = document.getElementById("rd-pipeline-track");
      if (track && !track.dataset.jumpBound) {
        track.dataset.jumpBound = "1";
        track.addEventListener("click", onJump);
      }
      const list = document.getElementById("rd-progress-list");
      if (list && !list.dataset.jumpBound) {
        list.dataset.jumpBound = "1";
        list.addEventListener("click", onJump);
      }
    },

    jumpToStep(id) {
      const targets = STEP_TARGETS[id] || ["#rd-results", ".rd-logs"];
      let dest = null;
      for (let i = 0; i < targets.length; i += 1) {
        const key = targets[i];
        const el = key.charAt(0) === "#" || key.charAt(0) === "."
          ? $(key)
          : document.getElementById(key);
        if (!el) {
          continue;
        }
        const shown = el.offsetParent !== null || getComputedStyle(el).display !== "none";
        if (shown && el.getClientRects().length) {
          dest = el;
          break;
        }
        if (!dest) {
          dest = el;
        }
      }
      if (dest) {
        this.scrollTo(dest);
      }
    },

    clearFailState() {
      this.failed = false;
      this.failedProc = "";
      this.failMessage = "";
      this.doneAtFail = null;
      this.stopping = false;
      this.frozenRows = null;
      this.frozenKey = null;
      this.nfFinished = false;
      this.endedAt = null;
      this.startedAt = null;
      this.launchAt = null;
    },

    /* —— Run lifecycle —— */
    bindRunHooks() {
      const runBt = document.getElementById("run-bt");
      if (runBt) {
        runBt.addEventListener("click", () => {
          this.clearFailState();
          this.startedAt = Date.now();
          this.launchAt = this.startedAt;
          this.stampRunId();
          // The form validates on submit, right after this click; a launched run disables the button.
          setTimeout(() => {
            if (!runBt.disabled && !this.failed) {
              // Click did not start a run (e.g. empty path caught before disable).
              this.launchAt = null;
              this.startedAt = null;
            }
            this.render();
          }, 0);
        });
      }
      const stopBt = document.getElementById("stop-bt");
      if (stopBt) {
        stopBt.addEventListener("click", () => {
          if (this.runState(this.procs()) === "running") {
            this.freeze();
          }
          this.stopping = true;
        });
      }
      const logs = document.querySelector(".rd-center .rd-logs");
      if (logs) {
        logs.style.height = "";
      }
    },

    hookTerminal() {
      const term = terminal();
      if (!term || term._rdHooked) {
        return;
      }
      term._rdHooked = true;
      const origAppend = term.append.bind(term);
      term.append = (raw, force) => {
        const text = String(raw || "");
        this.inspectLine(text);
        // Capture the detailed error body that accompanies "Failed to start".
        if (
          force &&
          /does not exist|No ONT FASTQ|No PacBio FASTQ|No FASTQ|FASTQ directory|Invalid FASTQ|Choose a FASTQ|Did you mean/i.test(
            text
          )
        ) {
          this.failMessage = cleanFailMessage(text);
          this.markFailed("");
        }
        return origAppend(raw, force);
      };
      const origStatus = term.setStatus.bind(term);
      term.setStatus = (text, kind) => {
        if (kind === "error") {
          const msg = cleanFailMessage(text);
          // Keep a detailed body message; do not replace it with the short status label.
          if (msg && !/Failed to start/i.test(msg) && !this.failMessage) {
            this.failMessage = msg;
          }
          if (/Failed to start/i.test(String(text || ""))) {
            this.markFailed("");
          }
        }
        return origStatus(text, kind);
      };
    },

    inspectLine(line) {
      if (/N E X T F L O W/i.test(line) || /Launching/i.test(line)) {
        const m = line.match(/version\s+([0-9][0-9.]*)/i);
        if (m) {
          this.setText("rd-nf-ver", m[1]);
        }
      }
      // Do not treat SSE "Process completed" / stream end as success — that raced
      // ahead of Nextflow and painted every stage green while assembly still ran.
      if (/BactFlow:\s*Nextflow finished successfully/i.test(line)) {
        this.markCompleted();
      }
      if (/Process failed|Nextflow exited with code\s+[1-9]|ERROR\s*~\s*\.nextflow\/history\.lock/i.test(line)) {
        this.failMessage = cleanFailMessage(line) || this.failMessage;
        this.markFailed("");
        return;
      }
      if (this.stopping || !FAIL_RE.test(line)) {
        return;
      }
      this.failMessage = cleanFailMessage(line) || this.failMessage;
      const m = line.match(/process\s*>\s*'?([A-Za-z_]\w*)/i);
      this.markFailed(m ? m[1] : "");
    },

    markFailed(procName) {
      if (!this.startedAt) {
        const first = this.procs().reduce((m, p) => Math.min(m, p.startedAt || Infinity), Infinity);
        this.startedAt = isFinite(first) ? first : Date.now();
      }
      if (!this.failed) {
        const eta = tracker();
        this.doneAtFail = new Set();
        if (eta) {
          eta.processes.forEach((p) => {
            if (p.done) {
              this.doneAtFail.add(p.id);
            }
          });
        }
      }
      this.failed = true;
      this.nfFinished = false;
      this.serverRunning = false;
      if (procName && !this.failedProc) {
        this.failedProc = norm(procName);
      }
      if (!this.endedAt) {
        this.endedAt = Date.now();
      }
      // Server-side: kill Nextflow/Java/tool orphans even if the UI already looks idle.
      if (typeof window.ensurePipelineStoppedOnFail === "function") {
        window.ensurePipelineStoppedOnFail();
      }
      if (typeof BactflowProcessEta !== "undefined" && typeof BactflowProcessEta.finalizeAll === "function") {
        BactflowProcessEta.finalizeAll(
          typeof BactflowTerminal !== "undefined" ? BactflowTerminal : null,
          "failed"
        );
      }
    },

    hookMeters() {
      const meters = window.BactflowMeters;
      if (!meters || meters._rdHooked) {
        return;
      }
      meters._rdHooked = true;
      const orig = meters.apply.bind(meters);
      meters.apply = (s) => {
        if (s && (s.cpu_percent != null || s.ram_percent != null)) {
          this.serverRunning = !!s.running;
          this.renderGauges(s);
        }
        return orig(s);
      };
      if (!meters.timer) {
        meters.tick();
        meters.timer = setInterval(() => meters.tick(), 1000);
      }
    },

    stampRunId() {
      const d = new Date();
      const pad = (n) => String(n).padStart(2, "0");
      this.setText(
        "rd-run-id",
        `#${d.getFullYear()}-${pad(d.getMonth() + 1)}-${pad(d.getDate())}_${pad(d.getHours())}-${pad(d.getMinutes())}`
      );
    },

    /* —— Planned stages from the form, in Nextflow execution order —— */
    val(id) {
      const el = document.getElementById(id);
      return el ? String(el.value || "").trim().toLowerCase() : "";
    },

    isTrue(id) {
      return this.val(id) === "true";
    },

    getPlannedStages() {
      const S = (id, name, icon, match, hint, side) => ({ id, name, icon, match, hint, side: !!side });
      const stages = [
        S("envsetup", "Environment setup", "setup", /^envsetup/, "envsetup"),
        S("testify", "Environment check", "check", /^testify/, "testify"),
      ];
      if (this.isTrue("setup_only")) {
        return stages;
      }

      const asm = this.val("assemblerDropdown");
      const concat = this.isTrue("concat_reads");
      if (asm === "pacbio") {
        if (concat) {
          stages.push(S("concat", "Read pooling", "pool", /^fastqconcater/, "concat_reads"));
        }
        stages.push(S("pbcheck", "PacBio read check", "readcheck", /^pacbioreadcheck/, "pacbio_read_check"));
        // Quality filter always before coverage when both are selected.
        if (this.isTrue("nanofilter")) {
          stages.push(S("readfilt", "Read filter", "filter", /^nanoreadfilt/, "readfilt"));
        }
        if (this.isTrue("coverage_filter")) {
          stages.push(S("covfilt", "Coverage filter", "coverage", /^coveragefilt/, "coverage_filt"));
        }
        stages.push(S("asm", "Flye (PacBio)", "pacbio", /^assemblypacbio/, "assembly_pacbio"));
      } else if (asm === "flye" || asm === "unicycler") {
        if (concat) {
          stages.push(S("concat", "Read pooling", "pool", /^fastqconcater/, "concat_reads"));
        }
        stages.push(S("dedup", "Deduplication", "dedup", /^deduper/, "deduper"));
        // Quality filter always before coverage when both are selected.
        if (this.isTrue("nanofilter")) {
          stages.push(S("readfilt", "Read filter", "filter", /^nanoreadfilt/, "readfilt"));
        }
        if (this.isTrue("coverage_filter")) {
          stages.push(S("covfilt", "Coverage filter", "coverage", /^coveragefilt/, "coverage_filt"));
        }
        stages.push(
          asm === "flye"
            ? S("asm", "Flye (ONT)", "flye", /^assemblyflye/, "assembly_flye1")
            : S("asm", "Unicycler", "unicycler", /^assemblyunicycler/, "assembly_unicycler")
        );
      } else if (asm === "spades") {
        stages.push(S("asm", "SPAdes", "spades", /^assemblyspades/, "assembly_spades"));
      }

      const hasInput =
        (asm && asm !== "none") || !!document.getElementById("postForm");
      if (!hasInput) {
        return stages;
      }
      if (this.isTrue("circle_genome")) {
        stages.push(S("circ", "Circulator", "circulator", /^circulator/, "circulator"));
      }
      if (this.isTrue("bakta_annot")) {
        stages.push(S("bakta", "Bakta", "bakta", /^(baktaannot|prokannot)/, "bakta", true));
      }
      if (this.isTrue("tax_class")) {
        stages.push(S("gtdb", "GTDB-Tk", "gtdb", /^taxonomygtdbtk/, "gtdbtk", true));
      }
      if (this.isTrue("run_checkm")) {
        stages.push(S("checkm", "CheckM", "checkm", /^checkmlineage/, "checkm", true));
      }
      if (this.isTrue("run_plasmids")) {
        stages.push(S("genomad", "geNomad", "genomad", /^plasmiddetect/, "genomad", true));
      }
      if (this.isTrue("run_quast")) {
        stages.push(S("quast", "QUAST", "quast", /^quastcheck/, "quast", true));
      }
      return stages;
    },

    procs() {
      const eta = tracker();
      if (!eta || !eta.processes) {
        return [];
      }
      return [...eta.processes.values()].filter((p) => p.id !== "_executor");
    },

    isDone(p) {
      if (!p.done) {
        return false;
      }
      return this.doneAtFail ? this.doneAtFail.has(p.id) : true;
    },

    /* Resolve every planned stage to a state, percent and remaining time. */
    resolveStages(planned, list) {
      const eta = tracker();
      return planned.map((stage) => {
        const procs = list.filter((p) => stage.match.test(norm(p.name)));
        const expected = eta ? eta.hintFor(stage.hint) : 180;
        const out = { stage, procs, state: "wait", pct: 0, remaining: expected, expected, elapsed: 0, cached: false };
        if (this.failed && this.failedProc && stage.match.test(this.failedProc)) {
          out.state = "fail";
          out.failDetail = this.failMessage || "";
        }
        if (!procs.length) {
          return out;
        }
        out.cached = procs.every((p) => p.cached);
        out.elapsed = Math.max(
          ...procs.map((p) => (p.done ? p.frozenElapsed || 0 : eta ? eta.elapsedSeconds(p) : 0))
        );
        if (out.state === "fail") {
          out.pct = Math.round(procs.reduce((a, p) => a + (eta ? eta.displayPercent(p) : 0), 0) / procs.length);
          out.remaining = 0;
          return out;
        }
        if (procs.every((p) => this.isDone(p))) {
          out.state = "done";
          out.pct = 100;
          out.remaining = 0;
          return out;
        }
        const active = procs.filter((p) => !this.isDone(p));
        const live = active.filter((p) => !(eta && eta.isWaiting(p)));
        if (this.failed) {
          out.state = "wait";
          out.remaining = 0;
          return out;
        }
        if (live.length) {
          out.state = "run";
          out.pct = Math.round(procs.reduce((a, p) => a + (this.isDone(p) ? 100 : eta.displayPercent(p)), 0) / procs.length);
          const rem = live.map((p) => eta.remainingSeconds(p)).filter((v) => v != null);
          out.remaining = rem.length ? Math.max(...rem) : null;
          out.overrun = rem.length === 0 && live.some((p) => {
            const e = eta.computeEta(p);
            return e && e.label === "overrun";
          });
        }
        return out;
      });
    },

    /* A stage has truly begun only when a matching process is running or done.
       Queued Nextflow rows (`[- ] process > circulator`) must NOT count — that
       used to advance the frontier and paint Flye green while it was still
       assembling. */
    stageBegun(r) {
      if (!r) {
        return false;
      }
      if (r.state === "run" || r.state === "done") {
        return true;
      }
      return (r.procs || []).some(
        (p) => p && !p.aborted && (p.done || p.cached || p.phase === "running")
      );
    },

    /* If a later stage has begun executing, earlier empty ones were skipped/cached. */
    fillSkipped(rows) {
      const eta = tracker();
      let lastStarted = -1;
      rows.forEach((r, i) => {
        if (!r.stage.side && this.stageBegun(r)) {
          lastStarted = i;
        }
      });
      rows.forEach((r, i) => {
        if (i < lastStarted && !r.stage.side && !r.procs.length && r.state === "wait") {
          r.state = "done";
          r.pct = 100;
          r.remaining = 0;
          r.skipped = true;
        }
      });
      // Later main-chain stage is executing ⇒ earlier main stages finished
      // (even if their ✔ line was missed). Never use merely-queued procs.
      let frontier = -1;
      rows.forEach((r, i) => {
        if (!r.stage.side && this.stageBegun(r)) {
          frontier = i;
        }
      });
      rows.forEach((r, i) => {
        if (frontier < 0 || i >= frontier || r.stage.side || r.state === "fail") {
          return;
        }
        // Protect a live stage unless a later main-chain stage is already
        // executing (DAG proves we moved on; ✔ was probably missed).
        const stillLive = (r.procs || []).some(
          (p) => p && !p.done && !p.aborted && p.phase === "running"
        );
        const laterRunning = rows.slice(i + 1).some(
          (r2) =>
            !r2.stage.side &&
            (r2.procs || []).some((p) => p && !p.done && !p.aborted && p.phase === "running")
        );
        if (stillLive && !laterRunning) {
          return;
        }
        r.procs.forEach((p) => {
          if (!p.done && eta && typeof eta.markDone === "function") {
            eta.markDone(p);
          }
        });
        if (r.state === "run" || r.state === "wait") {
          r.state = "done";
          r.pct = 100;
          r.remaining = 0;
          if (!r.procs.length) {
            r.skipped = true;
          }
        }
      });
      if (this.failed && !rows.some((r) => r.state === "fail")) {
        // Prefer a process that was mid-run; otherwise the next pending main step
        // (e.g. Deduplication when the FASTQ channel was empty / path wrong).
        let culprit = rows.find((r) => r.procs.length && r.state !== "done");
        if (!culprit) {
          const lastDone = rows.reduce((idx, r, i) => (r.state === "done" && !r.stage.side ? i : idx), -1);
          culprit =
            rows.find((r, i) => i > lastDone && !r.stage.side && r.state === "wait") ||
            rows.find((r) => !r.stage.side && r.state === "wait");
        }
        if (culprit) {
          culprit.state = "fail";
          culprit.failDetail = this.failMessage || "";
        }
      }
      rows.forEach((r) => {
        if (r.state === "fail" && this.failMessage && !r.failDetail) {
          r.failDetail = this.failMessage;
        }
      });
      return rows;
    },

    /* A step only runs once everything it depends on has finished; side branches do not block later steps. */
    applyGates(rows) {
      let chainDone = true;
      rows.forEach((r) => {
        // Never demote / rewind a stage Nextflow is actually executing.
        const live = (r.procs || []).some((p) => p && !p.done && !p.aborted && p.phase === "running");
        if (!chainDone && r.state === "run" && !live) {
          r.state = "wait";
          r.pct = 0;
          r.remaining = r.expected;
          r.procs.forEach((p) => {
            if (!p.done) {
              p.phase = "waiting";
              p.runningAt = null;
            }
          });
        }
        if (!r.stage.side) {
          chainDone = chainDone && r.state === "done";
        }
      });
      return rows;
    },

    freeze() {
      const rows = this.applyGates(this.fillSkipped(this.resolveStages(this.getPlannedStages(), this.procs())));
      this.frozenRows = rows.map((r) => ({
        ...r,
        state: r.state === "run" ? "stop" : r.state,
        remaining: 0,
      }));
      this.frozenKey = rows.map((r) => r.stage.id).join("|");
      this.endedAt = Date.now();
    },

    markCompleted() {
      if (this.failed || this.stopping) {
        return;
      }
      this.nfFinished = true;
      if (!this.endedAt) {
        this.endedAt = Date.now();
      }
      this.serverRunning = false;
      const eta = tracker();
      if (eta) {
        eta.processes.forEach((p) => {
          if (!p.done && !p.aborted && p.id !== "_executor" && typeof eta.markDone === "function") {
            eta.markDone(p);
          }
        });
        if (typeof eta.stopTick === "function") {
          eta.stopTick();
        }
      }
      if (window.BactflowTerminal && typeof BactflowTerminal.setStatus === "function") {
        BactflowTerminal.setStatus("All processes completed", "ok");
      }
    },

    runState(list) {
      if (this.failed) {
        return "failed";
      }
      if (this.stopping || this.frozenRows) {
        return "stopped";
      }
      if (this.nfFinished && !this.failed) {
        return "completed";
      }
      if (this.endedAt && !this.failed) {
        const unfinished = list.some((p) => !p.done && !p.aborted && p.id !== "_executor");
        if (!unfinished) {
          return "completed";
        }
      }
      const runBt = document.getElementById("run-bt");
      const localActive = list.some((p) => !p.done && !p.aborted) && !this.endedAt;
      const launching = this.launchAt && !this.endedAt && Date.now() - this.launchAt < 20000;
      if (this.serverRunning || localActive || launching || (runBt && runBt.disabled && this.startedAt && !this.endedAt)) {
        return "running";
      }
      if (this.endedAt && !this.failed) {
        return "completed";
      }
      return "idle";
    },

    /* —— Render —— */
    render() {
      const list = this.procs();
      const planned = this.getPlannedStages();
      const state = this.runState(list);
      if (state === "running" && !this.startedAt) {
        const first = list.reduce((m, p) => Math.min(m, p.startedAt || Infinity), Infinity);
        this.startedAt = isFinite(first) ? first : Date.now();
        this.endedAt = null;
      }
      if (state === "idle" && this.startedAt && !this.endedAt) {
        if (list.length) {
          this.endedAt = Date.now();
        } else {
          this.startedAt = null;
        }
      }
      let rows;
      if (this.frozenRows && this.frozenKey === planned.map((st) => st.id).join("|")) {
        rows = this.frozenRows;
      } else if (this.stopping || this.frozenRows) {
        if (!this.frozenRows) {
          this.freeze();
        }
        rows = this.frozenRows || this.applyGates(this.fillSkipped(this.resolveStages(planned, list)));
      } else {
        this.frozenRows = null;
        rows = this.applyGates(this.fillSkipped(this.resolveStages(planned, list)));
        // Only collapse remaining stages after a verified Nextflow success.
        // Never do this from a generic "completed" guess while serverRunning.
        if (this.nfFinished && !this.failed && !this.stopping && !this.serverRunning) {
          rows.forEach((r) => {
            if (r.state === "wait" || r.state === "run") {
              r.state = "done";
              r.pct = 100;
              r.remaining = 0;
              if (!r.procs.length) {
                r.skipped = true;
              }
            }
          });
        } else if (state === "running" && rows.length && rows.every((r) => !r.procs.length && r.state === "wait")) {
          rows[0].state = "run";
          rows[0].elapsed = this.startedAt ? (Date.now() - this.startedAt) / 1000 : 0;
        }
      }

      this.renderStatus(state);
      this.renderPipeline(rows);
      this.renderProgress(rows);
      this.renderOverview(rows, state);
      this.renderSummary(rows);
      this.renderTimes();
      this.renderGauges();
      this.maybeAutoscroll();
      if (
        (state === "completed" || state === "idle") &&
        this.endedAt &&
        !this.failed &&
        !this.stopping &&
        !this.frozenRows
      ) {
        if (typeof window.revealResultsOnSuccess === "function") {
          window.revealResultsOnSuccess(true);
        }
      }
    },

    alignLogsWithNextflow() {
      const logs = document.querySelector(".rd-center .rd-logs");
      const center = document.querySelector(".rd-center");
      const right = document.querySelector(".rd-right");
      const nf = right && right.querySelector(":scope > .rd-card:last-child");
      if (!logs || !center || !right || !nf || logs.classList.contains("is-hidden") || document.documentElement.classList.contains("rd-live-results")) {
        if (logs) {
          logs.style.height = "";
        }
        return;
      }
      if (window.matchMedia("(max-width: 1180px)").matches) {
        logs.style.height = "";
        return;
      }
      const cRect = center.getBoundingClientRect();
      const rRect = right.getBoundingClientRect();
      const nfBottom = nf.getBoundingClientRect().bottom - rRect.top + right.scrollTop;
      const logsTop = logs.getBoundingClientRect().top - cRect.top + center.scrollTop;
      const next = Math.max(220, Math.round(nfBottom - logsTop));
      const px = `${next}px`;
      if (logs.style.height !== px) {
        const body = document.getElementById("output-bactflow");
        const keep = body ? body.scrollTop : 0;
        const follow = !!(window.BactflowTerminal && window.BactflowTerminal.followScroll);
        logs.style.height = px;
        if (follow && window.BactflowTerminal && typeof window.BactflowTerminal.scrollToBottom === "function") {
          window.BactflowTerminal.scrollToBottom();
        } else if (body) {
          body.scrollTop = keep;
        }
      }
    },

    bindLogsAlign() {
      const onAlign = () => this.alignLogsWithNextflow();
      window.addEventListener("resize", onAlign);
      const right = document.querySelector(".rd-right");
      if (right && typeof ResizeObserver !== "undefined") {
        new ResizeObserver(onAlign).observe(right);
      }
      requestAnimationFrame(onAlign);
    },

    renderStatus(state) {
      const badge = document.getElementById("rd-run-badge");
      if (!badge) {
        return;
      }
      const label =
        state === "running"
          ? "Running"
          : state === "failed"
            ? "Failed"
            : state === "completed"
              ? "Done"
              : state === "stopped"
                ? "Stopped"
                : "Idle";
      const cls =
        state === "completed"
          ? "rd-status is-completed"
          : state === "stopped"
            ? "rd-status is-stopped"
            : "rd-status is-" + state;
      if (badge.className !== cls) {
        badge.className = cls;
      }
      if (badge.textContent !== label) {
        badge.textContent = label;
      }
    },

    renderPipeline(rows) {
      const track = document.getElementById("rd-pipeline-track");
      if (!track) {
        return;
      }
      const key = rows.map((r) => r.stage.id).join("|");
      if (key !== this.plannedKey) {
        this.plannedKey = key;
        track.innerHTML = "";
        rows.forEach((r) => {
          const el = document.createElement("div");
          el.className = "rd-stage is-wait";
          el.dataset.step = r.stage.id;
          el.innerHTML = '<div class="rd-stage-node"></div><div class="rd-stage-name"></div><div class="rd-stage-meta"></div>';
          el.querySelector(".rd-stage-name").textContent = r.stage.name;
          track.appendChild(el);
        });
      }
      const els = track.querySelectorAll(".rd-stage");
      rows.forEach((r, i) => {
        const el = els[i];
        if (!el) {
          return;
        }
        el.className = "rd-stage is-" + r.state;
        el.classList.toggle("is-jump", r.state === "done");
        el.title = r.state === "done" ? "Open " + r.stage.name : r.stage.name;
        const icon = r.state === "done" ? ICONS.done : r.state === "fail" ? ICONS.fail : ICONS[r.stage.icon];
        const node = el.querySelector(".rd-stage-node");
        if (node.dataset.icon !== icon) {
          node.innerHTML = icon;
          node.dataset.icon = icon;
        }
        let meta = "Pending";
        if (r.state === "done") {
          meta = r.skipped ? "Done" : r.cached ? "Cached" : "Done · " + this.formatDuration(r.elapsed);
        } else if (r.state === "fail") {
          meta = r.failDetail || this.failMessage || "Failed";
          el.title = meta;
        } else if (r.state === "stop") {
          meta = `Stopped at ${r.pct}%`;
        } else if (r.state === "run") {
          if (r.remaining != null && r.remaining > 0) {
            meta = `${r.pct}% · ${this.formatDuration(r.elapsed)} · ETA ${this.formatDuration(r.remaining)}`;
          } else if (r.overrun) {
            meta = `${r.pct}% · ${this.formatDuration(r.elapsed)} · still running`;
          } else {
            meta = `${r.pct}% · ${this.formatDuration(r.elapsed)} · running`;
          }
        } else if (r.procs.length) {
          meta = "Queued";
        }
        el.querySelector(".rd-stage-meta").textContent = meta;
      });
    },

    renderProgress(rows) {
      const box = document.getElementById("rd-progress-list");
      if (!box) {
        return;
      }
      if (box.childElementCount !== rows.length || box.dataset.key !== this.plannedKey) {
        box.dataset.key = this.plannedKey;
        box.innerHTML = "";
        rows.forEach((r) => {
          const row = document.createElement("div");
          row.dataset.step = r.stage.id;
          row.innerHTML = '<span class="rd-progress-dot"></span><span class="name"></span><span class="pct"></span><div class="track"><i></i></div>';
          row.querySelector(".name").textContent = r.stage.name;
          box.appendChild(row);
        });
      }
      rows.forEach((r, i) => {
        const row = box.children[i];
        row.className = "rd-progress-row is-" + r.state;
        row.classList.toggle("is-jump", r.state === "done");
        row.querySelector(".pct").textContent = r.state === "fail" ? "failed" : r.state === "stop" ? "stopped " + r.pct + "%" : r.pct + "%";
        row.querySelector(".track > i").style.width = r.pct + "%";
        if (r.state === "fail") {
          const detail = r.failDetail || this.failMessage || "Failed";
          row.title = detail;
          const nameEl = row.querySelector(".name");
          if (nameEl) {
            nameEl.textContent = r.stage.name;
          }
          let why = row.querySelector(".rd-fail-why");
          if (!why) {
            why = document.createElement("div");
            why.className = "rd-fail-why";
            row.appendChild(why);
          }
          why.textContent = detail;
        } else {
          const why = row.querySelector(".rd-fail-why");
          if (why) {
            why.remove();
          }
        }
      });
    },

    overallPercent(rows) {
      if (this.nfFinished && !this.failed && !this.stopping) {
        return 100;
      }
      const total = rows.reduce((a, r) => a + r.expected, 0);
      if (!total) {
        return 0;
      }
      const done = rows.reduce((a, r) => a + (r.expected * r.pct) / 100, 0);
      return Math.max(0, Math.min(100, Math.round((100 * done) / total)));
    },

    pipelineRemaining(rows, running) {
      let seq = 0;
      let par = 0;
      for (const r of rows) {
        if (r.state === "done" || r.state === "fail" || r.state === "stop") {
          continue;
        }
        if (running && r.state === "run" && (r.remaining == null || r.overrun)) {
          continue;
        }
        const rem = running && r.state === "run" ? Math.max(0, Number(r.remaining) || 0) : (r.expected || 0);
        if (r.stage && r.stage.side) {
          par = Math.max(par, rem);
        } else {
          seq += rem;
        }
      }
      const total = seq + par;
      if (running && rows.some((r) => r.state === "run") && total <= 0) {
        return null;
      }
      return total;
    },

    renderOverview(rows, state) {
      const pct = this.overallPercent(rows);
      const doneN = rows.filter((r) => r.state === "done").length;
      const runners = rows.filter((r) => r.state === "run").slice().sort((a, b) => (b.pct || 0) - (a.pct || 0));
      const current = runners[0] || rows.find((r) => r.state === "fail");
      const runCount = runners.length;
      const ring = document.getElementById("rd-task-ring");
      const bar = document.getElementById("rd-task-bar");
      if (ring) {
        ring.style.setProperty("--rd-pct", pct);
        ring.classList.toggle("is-failed", state === "failed");
        ring.classList.toggle("is-stopped", state === "stopped");
      }
      if (bar) {
        bar.style.width = pct + "%";
        bar.parentElement.classList.toggle("is-failed", state === "failed");
        bar.parentElement.classList.toggle("is-stopped", state === "stopped");
      }
      this.setText("rd-task-ring-pct", pct + "%");
      this.setText("rd-task-ring-sub", `${doneN} of ${rows.length} steps`);
      if (state === "running" && current && window.BactflowTerminal) {
        const extra = runCount > 1 ? ` (+${runCount - 1} parallel)` : "";
        BactflowTerminal.setStatus(`${current.stage.name} · ${current.pct}%${extra}`, "run");
      } else if (state === "completed" && window.BactflowTerminal) {
        BactflowTerminal.setStatus("All processes completed", "ok");
      }

      let title = "Ready";
      let desc = "Pipeline: " + rows.map((r) => r.stage.name).join(" → ");
      if (state === "running") {
        title = current ? current.stage.name : "Starting…";
        if (runCount > 1) {
          desc = `${runCount} steps running in parallel · showing ${current.stage.name} (${current.pct}%)`;
        } else {
          desc = current
            ? `Step ${rows.indexOf(current) + 1} of ${rows.length} is running`
            : "Waiting for Nextflow to schedule the first step";
        }
      } else if (state === "failed") {
        title = "Pipeline failed";
        desc =
          this.failMessage ||
          (current ? `${current.stage.name} failed. Check the live log for details.` : "Check the live log for details.");
      } else if (state === "stopped" || this.frozenRows) {
        title = "Stopped";
        desc = `Run stopped after ${doneN} of ${rows.length} steps.`;
      } else if (state === "completed" || (this.endedAt && rows.length && doneN === rows.length)) {
        title = "Completed";
        desc = `All ${rows.length} steps finished.`;
      }
      this.setText("rd-task-title", title);
      this.setText("rd-task-desc", desc);

      let remaining = null;
      if (state === "running" || (state === "idle" && !this.endedAt)) {
        remaining = this.pipelineRemaining(rows, state === "running");
      }
      const etaLabel =
        remaining == null || (state === "running" && remaining <= 0)
          ? state === "running" || !this.endedAt
            ? "—"
            : "0s"
          : this.formatDuration(remaining);
      this.setText("rd-m-eta", etaLabel);
      this.setText(
        "rd-m-eta-sub",
        remaining == null || remaining <= 0 || state !== "running"
          ? state === "running" && remaining == null
            ? "estimating…"
            : ""
          : "done ~" + this.formatClockAt(Date.now() + remaining * 1000)
      );
    },

    renderSummary(rows) {
      const count = (s) => rows.filter((r) => r.state === s).length;
      const runN = count("run");
      const doneN = count("done");
      const failN = count("fail");
      this.setText("rd-sum-total", String(rows.length));
      this.setText("rd-sum-done", String(doneN));
      this.setText("rd-sum-run", String(runN));
      this.setText("rd-sum-wait", String(rows.length - doneN - runN - failN));
      this.setText("rd-sum-fail", String(failN || (this.failed ? 1 : 0)));
    },

    renderTimes() {
      const started = this.startedAt ? this.formatClockAt(this.startedAt) : "—";
      const elapsed = this.startedAt
        ? this.formatDuration(((this.endedAt || Date.now()) - this.startedAt) / 1000)
        : "0s";
      this.setText("rd-header-started", started);
      this.setText("rd-header-elapsed", elapsed);
      this.setText("rd-m-started", started);
      this.setText("rd-m-elapsed", elapsed);
      this.setText(
        "rd-m-started-sub",
        this.startedAt ? new Date(this.startedAt).toLocaleDateString() : "not started yet"
      );
      this.setText(
        "rd-m-elapsed-sub",
        !this.startedAt ? "waiting for run" : this.endedAt ? "run finished" : "since start"
      );
    },

    requestedCores() {
      const launch =
        Number(window.__BACTFLOW_LAUNCH_CPUS__) ||
        Number((this.lastStats || {}).cpu_budget) ||
        0;
      const n = parseInt(((document.getElementById("cpus") || {}).value || "").trim(), 10);
      let asked = Number.isFinite(n) && n > 0 ? n : 0;
      if (!asked && launch > 0) {
        asked = launch;
      }
      if (launch > 0 && asked > 0) {
        return Math.min(asked, launch);
      }
      return asked || (launch > 0 ? launch : 0);
    },

    requestedMemoryGb() {
      const n = parseFloat(((document.getElementById("memory_gb") || {}).value || "").trim());
      return Number.isFinite(n) && n > 0 ? n : 0;
    },

    cpuBudgetLabel(_hostCount, budgetHint) {
      // Denominator is the session --cpus budget (like memory GB cap), not host cores.
      const budget =
        Number(window.__BACTFLOW_LAUNCH_CPUS__) ||
        Number(budgetHint) ||
        this.requestedCores() ||
        0;
      const asked = this.requestedCores() || budget;
      if (asked && budget) {
        return `${asked} of ${budget} cores`;
      }
      if (asked) {
        return `${asked} cores`;
      }
      if (budget) {
        return `${budget} cores`;
      }
      return "";
    },

    renderGauges(s) {
      if (s) {
        this.lastStats = s;
      } else {
        s = this.lastStats || {};
      }
      const asked = this.requestedCores();
      const jobCores = Number(s.job_cores);
      const jobLive = s.running && Number.isFinite(jobCores);
      let cpu = Number(s.cpu_percent) || 0;
      if (jobLive && asked) {
        cpu = (jobCores / asked) * 100;
      }
      cpu = Math.max(0, Math.min(100, cpu));
      const cap = this.requestedMemoryGb() || Number(s.ram_total_gb) || 0;
      const jobRam = Number(s.job_rss_gb);
      const usedRam = s.running && Number.isFinite(jobRam) ? jobRam : Number(s.ram_used_gb) || 0;
      const ram = cap > 0
        ? Math.max(0, Math.min(100, (usedRam / cap) * 100))
        : Math.max(0, Math.min(100, Number(s.ram_percent) || 0));
      const cpuRing = document.getElementById("rd-cpu-ring");
      const ramRing = document.getElementById("rd-ram-ring");
      if (cpuRing) {
        cpuRing.style.setProperty("--rd-pct", cpu);
      }
      if (ramRing) {
        ramRing.style.setProperty("--rd-pct", ram);
      }
      this.setText("rd-cpu-pct", Math.round(cpu) + "%");
      this.setText("rd-ram-pct", Math.round(ram) + "%");
      this.setText("rd-cpu-sub", this.cpuBudgetLabel(s.cpu_count, s.cpu_budget));
      this.setText("rd-ram-sub", cap ? `${Number(usedRam).toFixed(1)} / ${cap} GB` : "");
      this.setText("rd-m-res", `${Math.round(cpu)}% / ${Math.round(ram)}%`);
    },

    maybeAutoscroll() {
      if (window.BactflowTerminal && typeof window.BactflowTerminal.pinScroll === "function") {
        window.BactflowTerminal.pinScroll();
        return;
      }
      if (!this.autoScroll) {
        return;
      }
      const body = document.getElementById("output-bactflow");
      if (body) {
        const slack = 48;
        const near = (body.scrollHeight - body.scrollTop - body.clientHeight) <= slack;
        if (near) {
          body.scrollTop = body.scrollHeight;
        }
      }
    },

    /* —— Helpers —— */
    formatDuration(sec) {
      const s = Math.max(0, Math.round(Number(sec) || 0));
      const h = Math.floor(s / 3600);
      const m = Math.floor((s % 3600) / 60);
      const r = s % 60;
      if (h) {
        return `${h}h ${m}m`;
      }
      return m ? `${m}m ${r}s` : `${r}s`;
    },

    formatClockAt(ms) {
      return new Date(ms).toLocaleTimeString([], { hour: "2-digit", minute: "2-digit" });
    },

    setText(id, text) {
      const el = document.getElementById(id);
      if (el && el.textContent !== text) {
        el.textContent = text;
      }
    },
  };

  window.BactflowRunDashboard = RunDashboard;
  const boot = () => RunDashboard.init();
  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", boot);
  } else {
    boot();
  }
})();
