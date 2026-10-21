(() => {
  const metaRoot = document.querySelector('meta[name="bactflow-root"]');
  const scriptRoot = ((metaRoot && metaRoot.content) || "").replace(/\/$/, "");
  const path = window.location.pathname || "/";
  const inFrame = window.self !== window.top;

  const HUB_KEEP = [
    "/hub/",
    "/qc/",
    "/qc",
    "/assembly/",
    "/assembly",
    "/post/",
    "/post",
    "/wiki",
    "/run-lock",
    "/run-lock/reset",
    "/host-info",
  ];

  function keepAbsolute(url) {
    if (!url || typeof url !== "string") return true;
    if (!url.startsWith("/")) return true;
    if (url === "/") return true;
    if (scriptRoot && (url === scriptRoot || url.startsWith(scriptRoot + "/") || url.startsWith(scriptRoot + "?"))) {
      return true;
    }
    return HUB_KEEP.some((prefix) => url === prefix || url.startsWith(prefix + "/") || url.startsWith(prefix + "?"));
  }

  function rewritePath(url) {
    if (keepAbsolute(url)) return url;
    return scriptRoot + url;
  }

  function rewrite(url) {
    if (!scriptRoot || !url || typeof url !== "string") return url;
    try {
      if (url.startsWith("http://") || url.startsWith("https://")) {
        const u = new URL(url);
        if (u.origin !== window.location.origin) return url;
        return u.origin + rewritePath(u.pathname + u.search + u.hash);
      }
    } catch (err) {
      /* fall through */
    }
    return rewritePath(url);
  }

  if (scriptRoot) {
    const nativeFetch = window.fetch.bind(window);
    window.fetch = (input, init) => {
      if (typeof input === "string") return nativeFetch(rewrite(input), init);
      if (input instanceof Request) {
        try {
          const u = new URL(input.url, window.location.origin);
          if (u.origin === window.location.origin) {
            const pathOnly = u.pathname + u.search + u.hash;
            const next = rewrite(pathOnly);
            if (next !== pathOnly) return nativeFetch(new Request(next, input), init);
          }
        } catch (err) {
          /* leave the request untouched */
        }
      }
      return nativeFetch(input, init);
    };

    const XHR = window.XMLHttpRequest;
    if (XHR && XHR.prototype && XHR.prototype.open) {
      const nativeOpen = XHR.prototype.open;
      XHR.prototype.open = function (method, url, ...rest) {
        return nativeOpen.call(this, method, rewrite(url), ...rest);
      };
    }

    const NativeES = window.EventSource;
    if (NativeES) {
      function PatchedES(url, config) {
        return new NativeES(rewrite(url), config);
      }
      PatchedES.prototype = NativeES.prototype;
      PatchedES.CONNECTING = NativeES.CONNECTING;
      PatchedES.OPEN = NativeES.OPEN;
      PatchedES.CLOSED = NativeES.CLOSED;
      window.EventSource = PatchedES;
    }
  }

  const LAUNCH_IDS = ["run-bt", "install-bt", "illumina-filter-bt", "read-cat", "read-stat", "vis-bt"];

  function applyLockButtons(info, moduleId) {
    const busy = Boolean(info && info.running && info.module);
    const lockOthers = busy && info.module !== moduleId;
    const msg = busy
      ? `${info.label || "Another module"} is still running. Finish or stop it before launching another process.`
      : "";
    LAUNCH_IDS.forEach((id) => {
      const el = document.getElementById(id);
      if (!el) return;
      if (lockOthers) {
        el.classList.add("is-run-locked");
        el.setAttribute("aria-disabled", "true");
        el.title = msg;
      } else {
        el.classList.remove("is-run-locked");
        el.removeAttribute("aria-disabled");
        if (el.title === msg) el.removeAttribute("title");
      }
    });
    const mem = document.getElementById("memory_gb");
    if (mem && info) {
      const maxGb = Number(info.max_memory_gb || 0);
      const defGb = Number(info.default_memory_gb || 8);
      if (maxGb) mem.max = String(maxGb);
      const cur = Number(mem.value || 0);
      if (!cur) mem.value = String(Math.min(defGb, maxGb || defGb));
      else if (maxGb && cur > maxGb) mem.value = String(maxGb);
      const hint = document.getElementById("memory-hint");
      if (hint && maxGb) {
        const launch = Number(info.launch_memory_gb || 0);
        hint.textContent = launch
          ? `Hard cap for this run (launcher --memory ${launch}g). Host has ${info.host_mem_gb} GB; the form cannot exceed ${maxGb} GB.`
          : `Hard cap for this run. This machine has ${info.host_mem_gb} GB; values above ${maxGb} GB are reduced so the host is not exhausted.`;
      }
    }
    // Hard session cap from bactflow.sh --cpus. Form cannot exceed launch budget.
    const cpus = document.getElementById("cpus");
    if (cpus && info) {
      const hostC = Number(info.host_cpus || 0);
      const launchC = Number(info.launch_cpus || window.__BACTFLOW_LAUNCH_CPUS__ || 0);
      const defC = Number(info.default_cpus || launchC || 1);
      const maxC = launchC > 0 ? launchC : Number(info.max_cpus || hostC || defC || 0);
      if (maxC) {
        cpus.max = String(maxC);
        cpus.setAttribute("max", String(maxC));
      }
      if (launchC > 0) {
        window.__BACTFLOW_LAUNCH_CPUS__ = launchC;
        const cur = Number(cpus.value || 0);
        // Always clamp; default to the launcher budget on first apply.
        if (!cpus.dataset.budgetApplied || !cur || cur > launchC) {
          cpus.value = String(launchC);
        } else {
          cpus.value = String(Math.min(cur, launchC));
        }
        cpus.dataset.budgetApplied = "1";
        if (!cpus.dataset.capBound) {
          cpus.dataset.capBound = "1";
          cpus.addEventListener("input", () => {
            const lim = Number(window.__BACTFLOW_LAUNCH_CPUS__ || cpus.max || 0);
            const v = Number(cpus.value || 0);
            if (lim > 0 && v > lim) cpus.value = String(lim);
            if (v < 1) cpus.value = "1";
          });
        }
      } else if (!cpus.dataset.budgetApplied) {
        cpus.dataset.budgetApplied = "1";
        if (defC) cpus.value = String(Math.min(defC, maxC || defC));
      } else if (maxC && Number(cpus.value || 0) > maxC) {
        cpus.value = String(maxC);
      }
      const cpuHint = document.getElementById("cpus-hint");
      if (cpuHint) {
        cpuHint.textContent = launchC
          ? `Hard cap: --cpus ${launchC} (session budget ${launchC} cores${hostC ? `; host has ${hostC}` : ""}).`
          : hostC
            ? `This machine has ${hostC} cores.`
            : "";
      }
    }
  }

  // Module pages inside the hub iframe keep their fetch prefix, not a second top bar.
  if (inFrame) {
    const frameModule =
      path.startsWith("/qc") ? "qc" :
      path.startsWith("/assembly") ? "assembly" :
      path.startsWith("/post") ? "post" :
      "";
    window.addEventListener("message", (e) => {
      if (e.origin !== window.location.origin) return;
      if (e.data && e.data.type === "bactflow-lock") applyLockButtons(e.data.info, frameModule);
    });
    fetch("/run-lock", { cache: "no-store" })
      .then((res) => (res.ok ? res.json() : null))
      .then((info) => { if (info) applyLockButtons(info, frameModule); })
      .catch(() => {});
    setInterval(() => {
      fetch("/run-lock", { cache: "no-store" })
        .then((res) => (res.ok ? res.json() : null))
        .then((info) => { if (info) applyLockButtons(info, frameModule); })
        .catch(() => {});
    }, 2000);
    return;
  }

  const thisModule =
    path.startsWith("/qc") ? "qc" :
    path.startsWith("/assembly") ? "assembly" :
    path.startsWith("/post") ? "post" :
    "";

  const hubbed = Boolean(scriptRoot) || path === "/" || path.startsWith("/wiki") || path.startsWith("/qc") || path.startsWith("/assembly") || path.startsWith("/post");
  if (!hubbed) return;

  document.documentElement.classList.add("bf-hubbed");

  const VIEWS = [
    { href: "/#home", id: "home", label: "Home" },
    { href: "/#qc", id: "qc", label: "QC", src: "/qc/" },
    { href: "/#assembly", id: "assembly", label: "Assembly", src: "/assembly/" },
    { href: "/#post", id: "post", label: "Post-assembly", src: "/post/" },
    { href: "/#wiki", id: "wiki", label: "Wiki" },
    { href: "https://github.com/farhadm1990/bactflow", id: "github", label: "GitHub", external: true },
  ];

  function hashView() {
    const raw = (window.location.hash || "").replace(/^#/, "").toLowerCase();
    if (raw === "qc" || raw === "assembly" || raw === "post" || raw === "wiki") return raw;
    return "home";
  }

  function isHubShell() {
    return Boolean(document.body && document.body.classList.contains("bf-hub-shell"));
  }

  function showView(id) {
    if (!isHubShell()) {
      const target = id === "home" ? "/" : id === "wiki" ? "/#wiki" : "/" + id.replace(/^\/+/, "") + (id === "wiki" ? "" : "/");
      if (id === "wiki") window.location.href = "/#wiki";
      else if (id === "home") window.location.href = "/";
      else window.location.href = "/" + id + "/";
      return;
    }
    const valid = VIEWS.some((v) => v.id === id && !v.external) ? id : "home";
    document.querySelectorAll("[data-bf-view]").forEach((el) => {
      const on = el.getAttribute("data-bf-view") === valid;
      el.hidden = !on;
      if (on && el.tagName === "IFRAME" && !el.getAttribute("data-loaded")) {
        el.src = el.getAttribute("data-src") || "";
        el.setAttribute("data-loaded", "1");
      }
    });
    document.querySelectorAll(".bf-chrome-nav a").forEach((a) => {
      a.classList.toggle("is-active", a.dataset.view === valid);
    });
    const nextHash = valid === "home" ? "#home" : "#" + valid;
    if (window.location.hash !== nextHash) {
      history.replaceState(null, "", "/" + nextHash);
    }
  }

  function mountChrome() {
    if (!document.body) {
      document.addEventListener("DOMContentLoaded", mountChrome);
      return;
    }
    if (document.querySelector(".bf-chrome")) {
      if (isHubShell()) showView(hashView());
      return;
    }

    const bar = document.createElement("header");
    bar.className = "bf-chrome";
    bar.innerHTML = `
      <a class="bf-chrome-brand" href="/#home" data-view="home">
        <img src="/hub/static/images/logo.png" alt="">
        <div>
          <strong>BactFlow</strong>
          <span>From reads to insights</span>
        </div>
      </a>
      <nav class="bf-chrome-nav" aria-label="BactFlow modules"></nav>
      <div class="bf-chrome-aside">
        <div class="bf-chrome-status" id="bf-run-status">Ready</div>
        <button type="button" class="bf-chrome-reset" id="bf-run-reset" title="Clear stalled run locks and free the workspace">Reset</button>
      </div>
    `;
    const resetBtn = bar.querySelector("#bf-run-reset");
    if (resetBtn) {
      resetBtn.addEventListener("click", () => {
        resetWorkspace(resetBtn);
      });
    }
    const nav = bar.querySelector(".bf-chrome-nav");
    VIEWS.forEach((item) => {
      const a = document.createElement("a");
      a.href = item.href;
      a.textContent = item.label;
      if (item.external) {
        a.target = "_blank";
        a.rel = "noopener noreferrer";
        a.dataset.external = "1";
      } else {
        a.dataset.view = item.id;
      }
      nav.appendChild(a);
    });
    document.body.prepend(bar);

    document.addEventListener("click", (e) => {
      const a = e.target.closest("[data-view]");
      if (!a || !isHubShell()) return;
      const id = a.dataset.view;
      if (!id) return;
      e.preventDefault();
      showView(id);
    });
    window.addEventListener("hashchange", () => {
      if (isHubShell()) showView(hashView());
    });

    if (isHubShell()) showView(hashView());
    pollLock();
    setInterval(pollLock, 2000);
  }

  function applyLock(info) {
    const status = document.getElementById("bf-run-status");
    const resetBtn = document.getElementById("bf-run-reset");
    const busy = Boolean(info && info.running && info.module);
    const nextStatus = busy
      ? `${info.label || "A module"} is running`
      : "Ready";
    if (status) {
      const wasBusy = status.classList.contains("is-busy");
      if (wasBusy !== busy) {
        status.classList.toggle("is-busy", busy);
      }
      if (status.textContent !== nextStatus) {
        status.textContent = nextStatus;
      }
    }
    if (resetBtn) {
      const need = busy;
      if (resetBtn.classList.contains("is-needed") !== need) {
        resetBtn.classList.toggle("is-needed", need);
      }
      const nextTitle = busy
        ? `Clear stalled ${info.label || "job"} lock and free the workspace`
        : "Clear stalled run locks and free the workspace";
      if (resetBtn.title !== nextTitle) {
        resetBtn.title = nextTitle;
      }
    }
    applyLockButtons(info, thisModule);
    document.querySelectorAll("iframe[data-bf-view]").forEach((frame) => {
      try {
        frame.contentWindow.postMessage({ type: "bactflow-lock", info }, window.location.origin);
      } catch (err) {
        /* iframe may not be ready */
      }
    });
  }

  window.addEventListener("message", (e) => {
    if (e.origin !== window.location.origin) return;
    if (e.data && e.data.type === "bactflow-lock") applyLock(e.data.info);
  });

  async function pollLock() {
    try {
      const res = await fetch("/run-lock", { cache: "no-store" });
      if (!res.ok) return;
      applyLock(await res.json());
    } catch (err) {
      /* hub may be down when a module is launched standalone */
    }
  }

  function notifyModulesReset(info) {
    const payload = { type: "bactflow-workspace-reset", info: info || { running: false } };
    document.querySelectorAll("iframe[data-bf-view]").forEach((frame) => {
      try {
        frame.contentWindow.postMessage(payload, window.location.origin);
      } catch (err) {
        /* iframe may not be ready */
      }
    });
    // Also tell this window when a module page is open top-level (not in the hub shell).
    try {
      window.postMessage(payload, window.location.origin);
    } catch (err) {
      /* ignore */
    }
  }

  async function resetWorkspace(btn) {
    if (btn) {
      btn.disabled = true;
      btn.textContent = "Resetting…";
    }
    try {
      const res = await fetch("/run-lock/reset", { method: "POST", cache: "no-store" });
      const data = res.ok ? await res.json() : null;
      if (!res.ok) {
        return;
      }
      // Tell every module iframe to flip Run/Stop buttons and freeze the progress UI.
      notifyModulesReset(data);
      applyLock(data || { running: false });
    } catch (err) {
      /* silent — status poll will refresh the chrome badge */
    } finally {
      if (btn) {
        btn.disabled = false;
        btn.textContent = "Reset";
      }
      pollLock();
    }
  }

  mountChrome();
})();
