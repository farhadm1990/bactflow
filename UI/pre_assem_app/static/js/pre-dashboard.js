/* BactFlow Pre-assembly Dashboard — step states from the QC requests, live resources, theme */
(function () {
  const ICONS = {
    pool:
      '<svg viewBox="0 0 24 24"><path d="M4 6h16M4 12h16M4 18h16"/><path d="M8 3v6M16 9v6M11 15v6"/></svg>',
    trim:
      '<svg viewBox="0 0 24 24"><circle cx="6" cy="6" r="3"/><circle cx="6" cy="18" r="3"/><path d="M8.5 7.5L20 18M8.5 16.5L20 6"/></svg>',
    filter:
      '<svg viewBox="0 0 24 24"><path d="M4 5h16l-6 7v6l-4 2v-8z"/></svg>',
    cov:
      '<svg viewBox="0 0 24 24"><path d="M4 18V8M9 18V4M14 18v-7M19 18V10"/><path d="M3 20h18"/></svg>',
    stats:
      '<svg viewBox="0 0 24 24"><path d="M4 5h16v14H4z"/><path d="M4 10h16M10 10v9"/></svg>',
    plots:
      '<svg viewBox="0 0 24 24"><path d="M4 19V5M4 19h16M8 15v4M12 11v8M16 8v11"/></svg>',
    done: '<svg viewBox="0 0 24 24"><path d="M5 12l5 5L20 7"/></svg>',
    fail: '<svg viewBox="0 0 24 24"><path d="M6 6l12 12M18 6L6 18"/></svg>',
  };

  const STEPS = {
    pool: { name: "Initiate read cat", icon: "pool", desc: "Lists FASTQ files in the input directory" },
    concat: { name: "Concat reads", icon: "pool", desc: "Pools barcode or sample folder FASTQs into one file per sample" },
    trim: { name: "Size filter", icon: "trim", desc: "Removes FASTQ files below the size threshold", optional: true },
    ilfilter: { name: "Quality filter", icon: "filter", desc: "Filters Illumina reads by length and quality" },
    lrfilter: { name: "Read filter", icon: "filter", desc: "Filters ONT or PacBio reads by length and quality" },
    covfilt: { name: "Coverage filter", icon: "cov", desc: "Downsamples long reads to the target coverage with rasusa" },
    stats: { name: "Read statistics", icon: "stats", desc: "Read counts and lengths with SeqKit" },
    plots: { name: "Quality plots", icon: "plots", desc: "Read quality and length plots" },
    install: { name: "BactFlow install", icon: "pool", desc: "Creates the BactFlow conda environment" },
  };

  const STEP_TARGETS = {
    pool: ["fastq-ls-div", "rd-results"],
    concat: ["fastq-ls-div", "rd-results"],
    trim: ["fastq-ls-div", "rd-results"],
    ilfilter: ["vis-error", "rd-results"],
    lrfilter: ["vis-error", "rd-results"],
    covfilt: ["vis-error", "rd-results"],
    stats: ["seq-div", "rd-results"],
    plots: ["bar-qual-plot", "qual-plot", "rd-results"],
  };

  const PLOT_DIV_IDS = ["qual-plot", "qual-heat-plot", "qual-heat-pool-plot"];

  const ENDPOINTS = {
    "/ls-fastq": "pool",
    "/trim-list": "trim",
    "/filter-illumina": "ilfilter",
    "/filter-long-reads": "lrfilter",
    "/reads-stat": "stats",
    "/plot-qual": "plots",
    "/install-bactflow": "install",
  };

  /* A request can succeed at HTTP level and still report a failure in its JSON body. */
  function failureOf(id, ok, data) {
    if (!ok || data.error) {
      return data.error || data.message || "Request failed";
    }
    if (id === "pool" && !data.html_table) {
      return "No FASTQ files listed";
    }
    if (id === "trim" && data.status !== "completed") {
      return "Trimming did not complete";
    }
    if (id === "stats" && !data.html_output) {
      return "No statistics returned";
    }
    if (id === "plots" && !data.graph) {
      return "No plot data returned";
    }
    return "";
  }

  const PreDashboard = {
    steps: {},
    lastId: null,
    installed: null,
    plannedKey: null,
    tickTimer: null,

    init() {
      this.initTheme();
      this.bindUserMenu();
      this.hookMeters();
      this.hookPlotly();
      this.bindPlotResize();
      this.bindEmptyState();
      this.bindStickyPipeline();
      this.bindStepJumps();
      const platform = document.getElementById("read_platform");
      if (platform) {
        platform.addEventListener("change", () => this.resetSteps());
      }
      const concat = document.getElementById("concat_reads");
      if (concat) {
        concat.addEventListener("change", () => this.resetSteps());
      }
      const readFilter = document.getElementById("read_filter");
      if (readFilter) {
        readFilter.addEventListener("change", () => this.resetSteps());
      }
      const coverageFilter = document.getElementById("coverage_filter");
      if (coverageFilter) {
        coverageFilter.addEventListener("change", () => this.resetSteps());
      }
      this.render();
    },

    /* —— Theme —— */
    initTheme() {
      const btn = document.getElementById("rd-theme-bt");
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
        this.restylePlots();
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

    /* Process chain stays on screen and shrinks while the results column is scrolled. */
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
      };
      scroller.addEventListener("scroll", update, { passive: true });
      update();
      this.syncPipelineCompact = update;
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
      const targets = STEP_TARGETS[id] || ["rd-results"];
      let dest = targets[targets.length - 1];
      for (let i = 0; i < targets.length; i += 1) {
        const el = document.getElementById(targets[i]);
        if (!el) {
          continue;
        }
        const shown = el.offsetParent !== null || getComputedStyle(el).display !== "none";
        if (shown && el.getClientRects().length) {
          dest = targets[i];
          break;
        }
      }
      if (typeof window.jumpToSection === "function") {
        window.jumpToSection(dest);
        return;
      }
      const el = document.getElementById(dest);
      if (el) {
        el.scrollIntoView({ behavior: "smooth", block: "start" });
      }
    },

    bindPlotResize() {
      const pane = document.getElementById("pre-results-pane");
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
          try {
            Plotly.Plots.resize(el);
          } catch (e) { /* ignore */ }
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

    /* —— Step tracking from the requests each button makes —— */
    hookFetch() {
      const orig = window.fetch.bind(window);
      window.fetch = (input, init) => {
        let path = "";
        try {
          path = new URL(typeof input === "string" ? input : input.url, window.location.href).pathname;
        } catch (e) { /* ignore */ }
        const id = ENDPOINTS[path];
        const req = orig(input, init);
        if (path === "/bactflow-check") {
          req.then((res) => res.clone().json()).then((data) => this.onInstallCheck(data)).catch(() => {});
          return req;
        }
        if (!id) {
          return req;
        }
        // List FASTQ; when Concat Reads is on, pooling runs in the same request.
        if (path === "/ls-fastq") {
          const concatOn =
            this.formVal("read_platform") !== "illumina" &&
            this.formVal("concat_reads") === "true";
          this.start("pool");
          if (concatOn) {
            this.start("concat");
          }
          req
            .then((res) =>
              res
                .clone()
                .text()
                .then((text) => {
                  let data = {};
                  try {
                    data = text ? JSON.parse(text) : {};
                  } catch (e) {
                    return this.finishListFastq("Unexpected server response");
                  }
                  return this.finishListFastq(failureOf("pool", res.ok, data));
                })
            )
            .catch((err) => this.finishListFastq((err && err.message) || "Request failed"));
          return req;
        }
        // Long-read filter can run quality and/or coverage in one request.
        if (path === "/filter-long-reads") {
          const qualOn = this.formVal("read_filter") === "true";
          const covOn = this.formVal("coverage_filter") === "true";
          if (qualOn) {
            this.start("lrfilter");
          }
          if (covOn && !qualOn) {
            this.start("covfilt");
          } else if (covOn && qualOn) {
            this.steps.covfilt = {
              state: "wait",
              startedAt: null,
              endedAt: null,
              error: "",
              done: 0,
              total: 0,
            };
            this.render();
          }
          req
            .then((res) =>
              res
                .clone()
                .text()
                .then((text) => {
                  let data = {};
                  try {
                    data = text ? JSON.parse(text) : {};
                  } catch (e) {
                    return this.finishLongReadFilter("Unexpected server response");
                  }
                  return this.finishLongReadFilter(failureOf("lrfilter", res.ok, data));
                })
            )
            .catch((err) => this.finishLongReadFilter((err && err.message) || "Request failed"));
          return req;
        }
        this.start(id);
        req
          .then((res) =>
            res
              .clone()
              .text()
              .then((text) => {
                let data = {};
                try {
                  data = text ? JSON.parse(text) : {};
                } catch (e) {
                  return this.finish(id, "Unexpected server response");
                }
                if (id === "install") {
                  return res.ok && data.running ? null : this.finish(id, res.ok ? "" : failureOf(id, false, data));
                }
                return this.finish(id, failureOf(id, res.ok, data));
              })
          )
          .catch((err) => this.finish(id, (err && err.message) || "Request failed"));
        return req;
      };
    },

    finishListFastq(error) {
      const err = error || "";
      if (this.steps.pool && this.steps.pool.state === "run") {
        this.finish("pool", err);
      }
      if (this.steps.concat && this.steps.concat.state === "run") {
        this.finish("concat", err);
      }
    },

    finishLongReadFilter(error) {
      const err = error || "";
      if (this.steps.lrfilter && this.steps.lrfilter.state === "run") {
        this.finish("lrfilter", err);
      }
      const cov = this.steps.covfilt;
      if (!cov) {
        return;
      }
      if (err) {
        if (cov.state === "run") {
          this.finish("covfilt", err);
        } else if (cov.state === "wait") {
          cov.state = "fail";
          cov.error = err;
          cov.endedAt = Date.now();
          this.render();
        }
        return;
      }
      if (cov.state === "run") {
        this.finish("covfilt", "");
      } else if (cov.state === "wait") {
        cov.state = "done";
        cov.error = "";
        cov.endedAt = Date.now();
        cov.done = cov.total || 1;
        cov.total = cov.total || 1;
        this.render();
      }
    },

    onInstallCheck(data) {
      if (!data || data.status !== "success") {
        return;
      }
      this.installed = !!data.installed;
      if (this.installed && this.steps.install && this.steps.install.state === "run") {
        this.finish("install", "");
      }
      this.render();
    },

    start(id) {
      this.steps[id] = { state: "run", startedAt: Date.now(), endedAt: null, error: "", done: 0, total: 0 };
      this.lastId = id;
      this.ensureTick();
      this.pollProgress();
      this.render();
    },

    finish(id, error) {
      const st = this.steps[id];
      if (!st || st.state !== "run") {
        return;
      }
      if (id === "plots") {
        if (this.plotsReadyTimer) {
          clearInterval(this.plotsReadyTimer);
          this.plotsReadyTimer = null;
        }
        if (!error) {
          this.finishPlotsWhenReady();
          return;
        }
      }
      this.commitFinish(id, error);
    },

    plotsReady() {
      return PLOT_DIV_IDS.every((id) => {
        const el = document.getElementById(id);
        if (!el || el.style.display === "none") {
          return false;
        }
        return el.classList.contains("js-plotly-plot") && Array.isArray(el.data) && el.data.length > 0;
      });
    },

    finishPlotsWhenReady() {
      const st = this.steps.plots;
      if (!st || st.state !== "run") {
        return;
      }
      if (this.plotsReady()) {
        this.commitFinish("plots", "");
        return;
      }
      if (this.plotsReadyTimer) {
        return;
      }
      const started = Date.now();
      this.plotsReadyTimer = setInterval(() => {
        const cur = this.steps.plots;
        if (!cur || cur.state !== "run") {
          clearInterval(this.plotsReadyTimer);
          this.plotsReadyTimer = null;
          return;
        }
        if (this.plotsReady()) {
          clearInterval(this.plotsReadyTimer);
          this.plotsReadyTimer = null;
          this.commitFinish("plots", "");
        } else if (Date.now() - started > 180000) {
          clearInterval(this.plotsReadyTimer);
          this.plotsReadyTimer = null;
          this.commitFinish("plots", "");
        }
      }, 80);
    },

    commitFinish(id, error) {
      const st = this.steps[id];
      if (!st || st.state !== "run") {
        return;
      }
      st.state = error ? "fail" : "done";
      st.error = error || "";
      st.endedAt = Date.now();
      if (!error && st.total) {
        st.done = st.total;
      }
      this.render();
    },

    pollProgress() {
      if (this.progressTimer) {
        return;
      }
      const pull = () => {
        if (!this.activeIds().length) {
          if (this.progressTimer) {
            clearInterval(this.progressTimer);
            this.progressTimer = null;
          }
          return;
        }
        fetch("/pre_progress")
          .then((res) => res.json())
          .then((data) => {
            const job = data && data.job;
            if (job === "covfilt") {
              if (this.steps.lrfilter && this.steps.lrfilter.state === "run") {
                this.commitFinish("lrfilter", "");
              }
              if (!this.steps.covfilt || this.steps.covfilt.state === "wait") {
                this.start("covfilt");
              }
            }
            const st = job ? this.steps[job] : null;
            if (st && st.state === "run") {
              st.done = Number(data.done) || 0;
              st.total = Number(data.total) || 0;
              this.render();
            }
          })
          .catch(() => {});
      };
      pull();
      this.progressTimer = setInterval(pull, 400);
    },

    resetSteps() {
      Object.keys(this.steps).forEach((id) => {
        if (this.steps[id].state !== "run" && id !== "install") {
          delete this.steps[id];
        }
      });
      if (this.lastId && !this.steps[this.lastId]) {
        this.lastId = null;
      }
      this.render();
    },

    ensureTick() {
      if (!this.tickTimer) {
        this.tickTimer = setInterval(() => {
          this.render();
          if (!this.activeIds().length) {
            clearInterval(this.tickTimer);
            this.tickTimer = null;
          }
        }, 1000);
      }
    },

    activeIds() {
      return Object.keys(this.steps).filter((id) => this.steps[id].state === "run");
    },

    formVal(id) {
      return ((document.getElementById(id) || {}).value || "").toLowerCase();
    },

    plannedIds() {
      const platform = this.formVal("read_platform") || "auto";
      const illumina = platform === "illumina";
      const concatOn = !illumina && this.formVal("concat_reads") === "true";
      const filterOn = this.formVal("read_filter") === "true";
      const covOn = !illumina && this.formVal("coverage_filter") === "true";
      const ids = ["pool"];
      if (concatOn) {
        ids.push("concat");
      }
      if (filterOn) {
        ids.push(illumina ? "ilfilter" : "lrfilter");
      }
      if (covOn) {
        ids.push("covfilt");
      }
      if (!illumina && this.steps.trim) {
        ids.push("trim");
      }
      ids.push("stats", "plots");
      return ids;
    },

    rows() {
      return this.plannedIds().map((id) => {
        const st = this.steps[id] || { state: "wait", done: 0, total: 0 };
        const elapsed = st.startedAt ? ((st.endedAt || Date.now()) - st.startedAt) / 1000 : 0;
        const done = Number(st.done) || 0;
        const total = Number(st.total) || 0;
        let pct = 0;
        if (st.state === "done") {
          pct = 100;
        } else if (st.state === "run" && total > 0) {
          pct = Math.max(0, Math.min(99, Math.round((100 * done) / total)));
        } else if (total > 0) {
          pct = Math.max(0, Math.min(100, Math.round((100 * done) / total)));
        }
        return { id, def: STEPS[id], state: st.state, elapsed, error: st.error || "", done, total, pct };
      });
    },

    /* —— Render —— */
    render() {
      const rows = this.rows();
      this.renderPipeline(rows);
      this.renderProgress(rows);
      this.renderOverview(rows);
      this.renderSummary(rows);
    },

    runState() {
      if (this.activeIds().length) {
        return "running";
      }
      const last = this.lastId && this.steps[this.lastId];
      return last && last.state === "fail" ? "failed" : "idle";
    },

    renderPipeline(rows) {
      const track = document.getElementById("rd-pipeline-track");
      if (!track) {
        return;
      }
      const key = rows.map((r) => r.id).join("|");
      if (key !== this.plannedKey) {
        this.plannedKey = key;
        track.innerHTML = "";
        rows.forEach((r) => {
          const el = document.createElement("div");
          el.className = "rd-stage is-wait";
          el.dataset.step = r.id;
          el.innerHTML = '<div class="rd-stage-node"></div><div class="rd-stage-name"></div><div class="rd-stage-meta"></div>';
          el.querySelector(".rd-stage-name").textContent = r.def.name;
          el.title = r.def.desc;
          track.appendChild(el);
        });
      }
      const els = track.querySelectorAll(".rd-stage");
      rows.forEach((r, i) => {
        const el = els[i];
        if (!el) {
          return;
        }
        const complete = r.state === "done";
        el.className = "rd-stage is-" + (r.state === "fail" ? "fail" : complete ? "done" : r.state);
        el.classList.toggle("is-jump", complete && r.state !== "fail");
        el.title = complete && r.state !== "fail" ? "Open " + r.def.name : r.def.desc;
        const icon = complete && r.state !== "fail" ? ICONS.done : r.state === "fail" ? ICONS.fail : ICONS[r.def.icon];
        const node = el.querySelector(".rd-stage-node");
        if (node.dataset.icon !== icon) {
          node.innerHTML = icon;
          node.dataset.icon = icon;
        }
        let meta = r.def.optional ? "Optional" : "Pending";
        if (r.state === "run") {
          meta = r.id === "plots" && r.total && r.pct >= 99
            ? "Drawing plots…"
            : r.total
              ? `${r.done}/${r.total}`
              : "Running · " + this.formatDuration(r.elapsed);
        } else if (r.state === "done") {
          meta = "Done · " + this.formatDuration(r.elapsed);
        } else if (r.state === "fail") {
          meta = "Failed";
        }
        this.setText(el.querySelector(".rd-stage-meta"), meta);
      });
    },

    renderProgress(rows) {
      const box = document.getElementById("rd-progress-list");
      if (!box) {
        return;
      }
      if (box.dataset.key !== this.plannedKey) {
        box.dataset.key = this.plannedKey;
        box.innerHTML = "";
        rows.forEach((r) => {
          const row = document.createElement("div");
          row.dataset.step = r.id;
          row.innerHTML = '<span class="rd-progress-dot"></span><span class="name"></span><span class="pct"></span><div class="track"><i></i></div>';
          row.querySelector(".name").textContent = r.def.name;
          box.appendChild(row);
        });
      }
      rows.forEach((r, i) => {
        const row = box.children[i];
        const complete = r.state === "done";
        row.className = "rd-progress-row is-" + (r.state === "fail" ? "fail" : complete ? "done" : r.state);
        row.classList.toggle("is-jump", complete && r.state !== "fail");
        const label =
          r.state === "fail"
            ? "failed"
            : complete
              ? "100%"
              : r.id === "plots" && r.total && r.pct >= 99
                ? "drawing…"
                : r.total
                  ? `${r.done}/${r.total}`
                  : "0%";
        this.setText(row.querySelector(".pct"), label);
        row.querySelector(".track > i").style.width = (complete ? 100 : r.pct) + "%";
      });
    },

    renderOverview(rows) {
      const state = this.runState();
      const required = rows.filter((r) => !r.def.optional);
      const doneN = rows.filter((r) => r.state === "done").length;
      const pct = required.length
        ? Math.round(required.reduce((sum, r) => sum + r.pct, 0) / required.length)
        : 0;
      const ring = document.getElementById("rd-task-ring");
      const bar = document.getElementById("rd-task-bar");
      if (ring) {
        ring.style.setProperty("--rd-pct", pct);
        ring.classList.toggle("is-failed", state === "failed");
      }
      if (bar) {
        bar.style.width = pct + "%";
        bar.parentElement.classList.toggle("is-failed", state === "failed");
      }
      this.setText("rd-task-ring-pct", pct + "%");
      this.setText("rd-task-ring-sub", `${doneN} of ${rows.length} steps`);

      const lastId = this.lastId;
      const last = lastId ? this.steps[lastId] : null;
      const lastDef = lastId ? STEPS[lastId] : null;
      let title = "Ready";
      let desc = "Steps: " + rows.map((r) => r.def.name).join(" → ");
      if (state === "running" && lastDef) {
        title = lastDef.name;
        desc = lastDef.desc + "…";
      } else if (state === "failed" && lastDef) {
        title = lastDef.name + " failed";
        desc = last.error || "Check the message in the results below.";
      } else if (required.length && required.every((r) => r.state === "done")) {
        title = "All steps done";
        desc = "Reads are ready for assembly.";
      } else if (last && last.state === "done") {
        title = lastDef.name + " done";
        const next = rows.find((r) => r.state === "wait" && !r.def.optional);
        desc = next ? `Next: ${next.def.name}.` : "Results are shown below.";
      }
      this.setText("rd-task-title", title);
      this.setText("rd-task-desc", desc);

      const badge = document.getElementById("rd-run-badge");
      if (badge) {
        badge.className = "rd-status is-" + state;
        this.setText(badge, state === "running" ? "Running" : state === "failed" ? "Failed" : "Idle");
      }

      const started = last ? this.formatClockAt(last.startedAt) : "—";
      const elapsed = last ? this.formatDuration(((last.endedAt || Date.now()) - last.startedAt) / 1000) : "0s";
      this.setText("pre-header-action", lastDef ? lastDef.name : "—");
      this.setText("rd-header-started", started);
      this.setText("rd-header-elapsed", last ? elapsed : "—");
      this.setText("rd-m-started", started);
      this.setText("rd-m-started-sub", last ? new Date(last.startedAt).toLocaleDateString() : "not started yet");
      this.setText("rd-m-elapsed", elapsed);
      this.setText("rd-m-elapsed-sub", !last ? "no step run yet" : last.endedAt ? "last step" : "current step");
      this.setText("pre-m-last", lastDef ? lastDef.name : "—");
      this.setText(
        "pre-m-last-sub",
        !last ? "" : last.state === "run" ? "running" : last.state === "fail" ? "failed" : "completed"
      );
    },

    renderSummary(rows) {
      const count = (s) => rows.filter((r) => r.state === s).length;
      const doneN = count("done");
      const runN = count("run");
      const failN = count("fail");
      this.setText("rd-sum-total", String(rows.length));
      this.setText("rd-sum-done", String(doneN));
      this.setText("rd-sum-run", String(runN));
      this.setText("rd-sum-wait", String(rows.length - doneN - runN - failN));
      this.setText("rd-sum-fail", String(failN));
    },

    /* —— Resources —— */
    hookMeters() {
      const meters = window.BactflowMeters;
      if (!meters || meters._preHooked) {
        return;
      }
      meters._preHooked = true;
      const orig = meters.apply.bind(meters);
      meters.apply = (s) => {
        if (s && (s.cpu_percent != null || s.ram_percent != null)) {
          this.renderGauges(s);
        }
        return orig(s);
      };
    },

    requestedMemoryGb() {
      const n = parseFloat(((document.getElementById("memory_gb") || {}).value || "").trim());
      return Number.isFinite(n) && n > 0 ? n : 0;
    },

    renderGauges(s) {
      const cpu = Math.max(0, Math.min(100, Number(s.cpu_percent) || 0));
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
      this.setText("rd-cpu-sub", s.cpu_count ? `${s.cpu_count} cores` : "");
      this.setText("rd-ram-sub", cap ? `${Number(usedRam).toFixed(1)} / ${cap} GB` : "");
      this.setText("rd-m-res", `${Math.round(cpu)}% / ${Math.round(ram)}%`);
    },

    /* —— Plots follow the app theme on screen; downloads use a white background —— */
    plotColors(light) {
      return light
        ? { font: "#0f172a", grid: "rgba(15, 23, 42, 0.12)", paper: "rgba(0,0,0,0)" }
        : { font: "#e2e8f0", grid: "rgba(148, 163, 184, 0.18)", paper: "rgba(0,0,0,0)" };
    },

    themeLayout(layout, colors) {
      const out = { paper_bgcolor: colors.paper, plot_bgcolor: colors.paper, "font.color": colors.font };
      Object.keys(layout || {}).forEach((k) => {
        if (/^[xy]axis\d*$/.test(k)) {
          out[k + ".gridcolor"] = colors.grid;
          out[k + ".zerolinecolor"] = colors.grid;
          out[k + ".linecolor"] = colors.grid;
        }
      });
      return out;
    },

    isLight() {
      return document.documentElement.getAttribute("data-theme") === "light";
    },

    preparePlot(el) {
      return el && el.closest && el.closest("#pre-results-pane");
    },

    hookPlotly() {
      if (!window.Plotly || Plotly._preHooked) {
        return;
      }
      Plotly._preHooked = true;
      const origNew = Plotly.newPlot.bind(Plotly);
      Plotly.newPlot = (el, data, layout, config) => {
        const gd = typeof el === "string" ? document.getElementById(el) : el;
        return origNew(el, data, layout, config).then((plot) => {
          if (!this.preparePlot(gd)) {
            return plot;
          }
          return Plotly.relayout(gd, this.themeLayout(gd._fullLayout, this.plotColors(this.isLight()))).then(() => {
            this.finishPlotsWhenReady();
            return plot;
          });
        });
      };
      if (typeof window.downloadPlot === "function") {
        const origDownload = window.downloadPlot;
        window.downloadPlot = (plotId, filename, format) => {
          const gd = document.getElementById(plotId);
          if (!gd || !gd._fullLayout) {
            return origDownload(plotId, filename, format);
          }
          const white = this.themeLayout(gd._fullLayout, { font: "#1f2937", grid: "#e5e7eb", paper: "#ffffff" });
          Plotly.relayout(gd, white)
            .then(() => {
              origDownload(plotId, filename, format);
              return new Promise((resolve) => setTimeout(resolve, 800));
            })
            .finally(() => Plotly.relayout(gd, this.themeLayout(gd._fullLayout, this.plotColors(this.isLight()))));
        };
      }
    },

    restylePlots() {
      if (!window.Plotly) {
        return;
      }
      const colors = this.plotColors(this.isLight());
      document.querySelectorAll("#pre-results-pane .js-plotly-plot").forEach((gd) => {
        Plotly.relayout(gd, this.themeLayout(gd._fullLayout, colors));
      });
    },

    /* Shows the hint only while no result block is visible. */
    bindEmptyState() {
      const pane = document.getElementById("pre-results-pane");
      const hint = document.getElementById("pre-results-empty");
      if (!pane || !hint) {
        return;
      }
      const update = () => {
        const any = [...pane.children].some((c) => c.style.display !== "none" && getComputedStyle(c).display !== "none");
        hint.style.display = any ? "none" : "";
      };
      new MutationObserver(update).observe(pane, { subtree: true, attributes: true, attributeFilter: ["style"] });
      update();
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

    setText(target, text) {
      const el = typeof target === "string" ? document.getElementById(target) : target;
      if (el && el.textContent !== text) {
        el.textContent = text;
      }
    },
  };

  window.BactflowPreDashboard = PreDashboard;
  // The install check runs on DOMContentLoaded from pre-assem.js, before init(), so watch requests right away.
  PreDashboard.hookFetch();
  const boot = () => PreDashboard.init();
  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", boot);
  } else {
    boot();
  }
})();
