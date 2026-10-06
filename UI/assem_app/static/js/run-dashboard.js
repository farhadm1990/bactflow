/* BactFlow Run Dashboard — nav + live sync from terminal/process ETA */
(function () {
  const $ = (sel, root) => (root || document).querySelector(sel);
  const $$ = (sel, root) => Array.from((root || document).querySelectorAll(sel));

  const RunDashboard = {
    startedAt: null,
    autoScroll: true,

    init() {
      this.bindNav();
      this.bindLogTools();
      this.bindRunHooks();
      this.stampRunId();
      this.setStatus("ready");
      setInterval(() => this.tick(), 1000);
      this.tick();
      const out = document.getElementById("output-div");
      if (out) {
        out.style.display = "flex";
      }
    },

    bindRunHooks() {
      const runBt = document.getElementById("run-bt");
      if (runBt) {
        runBt.addEventListener("click", () => {
          this.startedAt = Date.now();
          this.stampRunId();
          this.setStatus("running");
          this.showView("overview");
        });
      }
      const repBt = document.getElementById("rep-bt");
      if (repBt) {
        repBt.addEventListener("click", () => {
          setTimeout(() => this.showView("results"), 50);
        });
      }
    },

    bindNav() {
      $$(".rd-nav-btn[data-view]").forEach((btn) => {
        btn.addEventListener("click", () => this.showView(btn.dataset.view));
      });
      const viewLogs = $("#rd-view-logs-btn");
      if (viewLogs) {
        viewLogs.addEventListener("click", () => {
          this.showView("overview");
          const logs = $(".rd-logs");
          if (logs) {
            logs.scrollIntoView({ behavior: "smooth", block: "nearest" });
          }
        });
      }
    },

    showView(name) {
      $$(".rd-nav-btn[data-view]").forEach((b) => {
        b.classList.toggle("is-active", b.dataset.view === name);
      });
      $$(".rd-view").forEach((v) => {
        v.classList.toggle("is-active", v.id === "rd-view-" + name);
      });
    },

    bindLogTools() {
      const auto = $("#rd-autoscroll");
      if (auto) {
        auto.checked = true;
        auto.addEventListener("change", () => {
          this.autoScroll = !!auto.checked;
        });
      }
      const dl = $("#rd-download-log");
      if (dl) {
        dl.addEventListener("click", () => {
          const body = document.getElementById("output-bactflow");
          const text = body ? body.innerText || body.textContent || "" : "";
          const blob = new Blob([text], { type: "text/plain" });
          const a = document.createElement("a");
          a.href = URL.createObjectURL(blob);
          a.download = "bactflow-log.txt";
          a.click();
          URL.revokeObjectURL(a.href);
        });
      }
      $$(".rd-tab[data-logtab]").forEach((tab) => {
        tab.addEventListener("click", () => {
          $$(".rd-tab[data-logtab]").forEach((t) => t.classList.toggle("is-active", t === tab));
        });
      });
    },

    stampRunId() {
      const el = $("#rd-run-id");
      if (!el) {
        return;
      }
      const d = new Date();
      const pad = (n) => String(n).padStart(2, "0");
      el.textContent =
        "#" +
        d.getFullYear() +
        "-" +
        pad(d.getMonth() + 1) +
        "-" +
        pad(d.getDate()) +
        "_" +
        pad(d.getHours()) +
        "-" +
        pad(d.getMinutes());
    },

    setStatus(kind) {
      const badge = $("#rd-run-badge");
      if (!badge) {
        return;
      }
      badge.className = "rd-badge is-" + (kind || "ready");
      badge.textContent =
        kind === "running" ? "Running" : kind === "failed" ? "Failed" : kind === "done" ? "Completed" : "Ready";
      const pulse = $("#rd-gpu-label");
      if (pulse && kind === "running") {
        pulse.textContent = "Pipeline running · Docker + Nextflow";
      }
    },

    formatDuration(sec) {
      const s = Math.max(0, Math.round(Number(sec) || 0));
      const m = Math.floor(s / 60);
      const r = s % 60;
      if (m >= 60) {
        const h = Math.floor(m / 60);
        return `${h}h ${m % 60}m`;
      }
      return m ? `${m}m ${r}s` : `${r}s`;
    },

    formatClock(secFromNow) {
      if (secFromNow == null || !isFinite(secFromNow)) {
        return "—";
      }
      const d = new Date(Date.now() + secFromNow * 1000);
      return d.toLocaleTimeString([], { hour: "2-digit", minute: "2-digit" }) + " CEST";
    },

    getEta() {
      const eta = window.BactflowProcessEta;
      if (!eta || !eta.processes) {
        return null;
      }
      return eta;
    },

    procs() {
      const eta = this.getEta();
      if (!eta) {
        return [];
      }
      return [...eta.processes.values()].filter((p) => p.id !== "_executor");
    },

    activeProc(list) {
      return (
        list.find((p) => !p.done && !etaWaiting(p, this.getEta())) ||
        list.find((p) => !p.done) ||
        list[list.length - 1] ||
        null
      );
    },

    tick() {
      const eta = this.getEta();
      const list = this.procs();
      const running = list.some((p) => !p.done);
      const any = list.length > 0;

      if (running) {
        if (!this.startedAt) {
          this.startedAt = Date.now();
        }
        this.setStatus("running");
      } else if (any && list.every((p) => p.done)) {
        this.setStatus("done");
      }

      this.renderPipeline(list, eta);
      this.renderTask(list, eta);
      this.renderSummary(list, eta);
      this.renderElapsed();
      this.syncGauges();
      this.detectNextflow();
      this.maybeAutoscroll();
    },

    detectNextflow() {
      const el = document.getElementById("rd-nf-ver");
      const body = document.getElementById("output-bactflow");
      if (!el || !body) {
        return;
      }
      const text = body.innerText || body.textContent || "";
      const m = text.match(/N E X T F L O W\s*[~\s]*version\s*([0-9.]+)/i) || text.match(/Nextflow\s+([0-9]+\.[0-9.]+)/i);
      if (m) {
        el.textContent = m[1];
      }
    },

    renderPipeline(list, eta) {
      const stages = $$(".rd-stage");
      if (!stages.length) {
        return;
      }
      // Map chronological process list onto stage slots; keep extras collapsed into last active.
      stages.forEach((stage, idx) => {
        const proc = list[idx];
        const nameEl = stage.querySelector(".rd-stage-name");
        const metaEl = stage.querySelector(".rd-stage-meta");
        const nodeEl = stage.querySelector(".rd-stage-node");
        if (!proc) {
          stage.className = "rd-stage is-wait";
          if (nameEl && !stage.dataset.fixed) {
            nameEl.textContent = idx === 0 ? "Waiting" : "—";
          }
          if (metaEl) {
            metaEl.textContent = "Waiting";
          }
          if (nodeEl) {
            nodeEl.textContent = "○";
          }
          return;
        }
        const waiting = etaWaiting(proc, eta);
        const done = !!proc.done;
        const pct = eta && eta.displayPercent ? eta.displayPercent(proc) : 0;
        stage.className = "rd-stage " + (done ? "is-done" : waiting ? "is-wait" : "is-run");
        if (nameEl) {
          nameEl.textContent = eta && eta.prettyName ? eta.prettyName(proc.name) : proc.name;
        }
        if (nodeEl) {
          nodeEl.textContent = done ? "✓" : waiting ? "○" : "◎";
        }
        if (metaEl) {
          if (done) {
            const sec = proc.frozenElapsed || 0;
            metaEl.textContent = `Completed ${this.formatDuration(sec)}${proc.cached ? " · cached" : ""}`;
          } else if (waiting) {
            metaEl.textContent = "Waiting";
          } else {
            const rem = eta && eta.remainingSeconds ? eta.remainingSeconds(proc) : null;
            metaEl.textContent =
              `${pct}%` + (rem != null ? ` · ~${this.formatDuration(rem)}` : "");
          }
        }
      });
    },

    renderTask(list, eta) {
      const active = this.activeProc(list);
      const title = $("#rd-task-title");
      const desc = $("#rd-task-desc");
      const ring = $("#rd-task-ring");
      const ringPct = $("#rd-task-ring-pct");
      const bar = $("#rd-task-bar");
      const listEl = $("#rd-task-list");
      const depsEl = $("#rd-task-deps");

      if (!active) {
        if (title) {
          title.textContent = "Ready";
        }
        if (desc) {
          desc.textContent = "Start a run from Analysis Setup to populate live progress.";
        }
        if (ring) {
          ring.style.setProperty("--rd-pct", 0);
        }
        if (ringPct) {
          ringPct.textContent = "0%";
        }
        if (bar) {
          bar.style.width = "0%";
        }
        this.setMetric("rd-m-elapsed", "—");
        this.setMetric("rd-m-eta", "—");
        this.setMetric("rd-m-total", "0");
        this.setMetric("rd-m-done", "0");
        this.setMetric("rd-m-run", "0");
        this.setMetric("rd-m-fail", "0");
        if (listEl) {
          listEl.innerHTML = "";
        }
        if (depsEl) {
          depsEl.innerHTML = "";
        }
        return;
      }

      const pct = eta && eta.displayPercent ? eta.displayPercent(active) : 0;
      const rem = eta && eta.remainingSeconds ? eta.remainingSeconds(active) : null;
      const elapsed = eta && eta.elapsedSeconds ? eta.elapsedSeconds(active) : 0;
      const pretty = eta && eta.prettyName ? eta.prettyName(active.name) : active.name;

      if (title) {
        title.textContent = pretty + (active.tag ? ` (${active.tag})` : "");
      }
      if (desc) {
        desc.textContent = active.done
          ? "Process finished"
          : etaWaiting(active, eta)
            ? "Queued — waiting for executor"
            : "Assembly / analysis process in progress";
      }
      if (ring) {
        ring.style.setProperty("--rd-pct", pct);
      }
      if (ringPct) {
        ringPct.textContent = pct + "%";
      }
      if (bar) {
        bar.style.width = pct + "%";
      }

      const doneN = list.filter((p) => p.done).length;
      const runN = list.filter((p) => !p.done && !etaWaiting(p, eta)).length;
      const waitN = list.filter((p) => !p.done && etaWaiting(p, eta)).length;

      this.setMetric("rd-m-elapsed", this.formatDuration(elapsed));
      this.setMetric("rd-m-eta", rem != null ? this.formatDuration(rem) : "—");
      this.setMetric("rd-m-total", String(list.length));
      this.setMetric("rd-m-done", String(doneN));
      this.setMetric("rd-m-run", String(runN));
      this.setMetric("rd-m-fail", "0");

      if (listEl) {
        listEl.innerHTML = "";
        list.forEach((p) => {
          const pPct = eta && eta.displayPercent ? eta.displayPercent(p) : 0;
          const waiting = etaWaiting(p, eta);
          const row = document.createElement("div");
          row.className = "rd-task-row " + (p.done ? "is-done" : waiting ? "is-wait" : "is-run");
          const label = (eta && eta.prettyName ? eta.prettyName(p.name) : p.name) + (p.tag ? ` (${p.tag})` : "");
          const status = p.done
            ? `Done · ${this.formatDuration(p.frozenElapsed || 0)}`
            : waiting
              ? "Waiting"
              : `Running · ${this.formatDuration(eta.elapsedSeconds(p))} / ${
                  rem != null && p === active ? this.formatDuration(rem) : "…"
                }`;
          row.innerHTML = `<span class="label"></span><span class="meta"></span><div class="rd-mini-bar"><i style="width:${p.done ? 100 : pPct}%"></i></div>`;
          row.querySelector(".label").textContent = label;
          row.querySelector(".meta").textContent = status;
          listEl.appendChild(row);
        });
      }

      if (depsEl) {
        const idx = list.indexOf(active);
        const deps = list.slice(idx + 1).filter((p) => !p.done);
        depsEl.innerHTML = "";
        deps.slice(0, 4).forEach((p) => {
          const row = document.createElement("div");
          row.className = "rd-task-row is-wait";
          row.innerHTML = `<span class="label"></span><span class="meta"></span><div class="rd-mini-bar"><i style="width:0%"></i></div>`;
          row.querySelector(".label").textContent = eta.prettyName ? eta.prettyName(p.name) : p.name;
          row.querySelector(".meta").textContent = "Waiting for " + pretty;
          depsEl.appendChild(row);
        });
      }

      const etaBig = $("#rd-eta-big");
      const etaSub = $("#rd-eta-sub");
      if (etaBig) {
        etaBig.textContent = rem != null ? "~ " + this.formatDuration(rem) : "—";
      }
      if (etaSub) {
        etaSub.textContent = rem != null ? this.formatClock(rem) : "Estimated completion";
      }
    },

    renderSummary(list) {
      const doneN = list.filter((p) => p.done).length;
      const runN = list.filter((p) => !p.done && !etaWaiting(p, this.getEta())).length;
      const waitN = list.filter((p) => !p.done && etaWaiting(p, this.getEta())).length;
      const total = Math.max(list.length, 1);
      this.setText("rd-sum-done", String(doneN));
      this.setText("rd-sum-run", String(runN));
      this.setText("rd-sum-wait", String(waitN));
      this.setText("rd-sum-fail", "0");
      const setW = (id, n) => {
        const el = $(id);
        if (el) {
          el.style.width = Math.round((100 * n) / total) + "%";
        }
      };
      setW("#rd-sum-done-bar", doneN);
      setW("#rd-sum-run-bar", runN);
      setW("#rd-sum-wait-bar", waitN);
    },

    renderElapsed() {
      const el = $("#rd-header-elapsed");
      if (!el) {
        return;
      }
      if (!this.startedAt) {
        el.textContent = "Elapsed —";
        return;
      }
      const sec = (Date.now() - this.startedAt) / 1000;
      el.textContent = "Elapsed " + this.formatDuration(sec);
    },

    syncGauges() {
      const cpuVal = document.getElementById("bf-cpu-val");
      const ramVal = document.getElementById("bf-ram-val");
      const cpuPct = parseFloat((cpuVal && cpuVal.textContent) || "") || 0;
      let ramPct = 0;
      if (ramVal && ramVal.textContent) {
        const m = ramVal.textContent.match(/\((\d+(?:\.\d+)?)%\)/);
        ramPct = m ? parseFloat(m[1]) : parseFloat(ramVal.textContent) || 0;
      }
      const cpuRing = $("#rd-cpu-ring");
      const ramRing = $("#rd-ram-ring");
      const cpuLabel = $("#rd-cpu-pct");
      const ramLabel = $("#rd-ram-pct");
      if (cpuRing) {
        cpuRing.style.setProperty("--rd-pct", Math.min(100, cpuPct));
      }
      if (ramRing) {
        ramRing.style.setProperty("--rd-pct", Math.min(100, ramPct));
      }
      if (cpuLabel) {
        cpuLabel.textContent = Math.round(cpuPct) + "%";
      }
      if (ramLabel) {
        ramLabel.textContent = Math.round(ramPct) + "%";
      }
    },

    maybeAutoscroll() {
      if (!this.autoScroll) {
        return;
      }
      const body = document.getElementById("output-bactflow");
      if (body) {
        body.scrollTop = body.scrollHeight;
      }
    },

    setMetric(id, text) {
      const el = document.getElementById(id);
      if (el) {
        el.textContent = text;
      }
    },

    setText(id, text) {
      const el = document.getElementById(id);
      if (el) {
        el.textContent = text;
      }
    },
  };

  function etaWaiting(proc, eta) {
    if (!proc || proc.done) {
      return false;
    }
    if (eta && typeof eta.isWaiting === "function") {
      return eta.isWaiting(proc);
    }
    return false;
  }

  window.BactflowRunDashboard = RunDashboard;
  document.addEventListener("DOMContentLoaded", () => RunDashboard.init());
})();
