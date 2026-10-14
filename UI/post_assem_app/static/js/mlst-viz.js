/* MLST table + Plotly ST / allele charts for Post-assembly */
(function () {
  const CHART_HEIGHT = 320;
  const PNG_RES_KEY = "bactflow_mlst_png_res";
  const PNG_DIMS_KEY = "bactflow_mlst_png_dims";
  const PNG_LEVELS = [
    { id: "basic", label: "Basic", scale: 1, cap: 2048 },
    { id: "medium", label: "Medium", scale: 2, cap: 4096 },
    { id: "high", label: "High", scale: 8 / 3, cap: 5461 },
    { id: "higher", label: "Higher", scale: 10 / 3, cap: 6827 },
    { id: "max", label: "Max", scale: 4, cap: 8192 },
  ];

  function plotColors() {
    const light = document.documentElement.getAttribute("data-theme") === "light";
    return light
      ? { font: "#1f2937", grid: "#e5e7eb", paper: "#ffffff", bar: "#7c3aed" }
      : { font: "#e5e7eb", grid: "#334155", paper: "rgba(0,0,0,0)", bar: "#a78bfa" };
  }

  function purge(el) {
    if (!el || !window.Plotly) {
      return;
    }
    try {
      Plotly.purge(el);
    } catch (e) {
      /* ignore */
    }
    el.innerHTML = "";
    el.style.height = CHART_HEIGHT + "px";
    el.style.maxHeight = CHART_HEIGHT + "px";
    el.style.overflow = "hidden";
  }

  function parseDim(value, fallback) {
    const n = Number(value);
    if (!Number.isFinite(n) || n <= 0) {
      return fallback;
    }
    return Math.max(240, Math.min(8192, Math.round(n)));
  }

  function pngDims(el) {
    const wEl = document.getElementById("mlst-png-w");
    const hEl = document.getElementById("mlst-png-h");
    const fallbackW = Math.max(640, (el && el.clientWidth) || 960);
    const fallbackH = Math.max(
      280,
      (el && el._fullLayout && el._fullLayout.height) || CHART_HEIGHT
    );
    return {
      width: parseDim(wEl && wEl.value, fallbackW),
      height: parseDim(hEl && hEl.value, fallbackH),
    };
  }

  function pngLevel() {
    const el = document.getElementById("mlst-png-res");
    const idx = el
      ? Math.max(0, Math.min(PNG_LEVELS.length - 1, Number(el.value) || 0))
      : 1;
    return PNG_LEVELS[idx] || PNG_LEVELS[1];
  }

  function downloadChart(elId, filename) {
    const el = document.getElementById(elId);
    if (!el || !window.Plotly || !el.data) {
      return Promise.resolve();
    }
    const dims = pngDims(el);
    const level = pngLevel();
    const scale = Math.max(
      1,
      Math.min(level.scale, level.cap / Math.max(1, dims.width), level.cap / Math.max(1, dims.height))
    );
    return Plotly.downloadImage(el, {
      format: "png",
      filename: filename.replace(/\.png$/i, ""),
      width: dims.width,
      height: dims.height,
      scale,
    });
  }

  function bindDownloadControls() {
    const slider = document.getElementById("mlst-png-res");
    const label = document.getElementById("mlst-png-res-label");
    const wEl = document.getElementById("mlst-png-w");
    const hEl = document.getElementById("mlst-png-h");

    if (slider && !slider.dataset.bound) {
      slider.dataset.bound = "1";
      const apply = (store) => {
        const idx = Math.max(0, Math.min(PNG_LEVELS.length - 1, Number(slider.value) || 0));
        slider.value = String(idx);
        const level = PNG_LEVELS[idx];
        if (label) {
          label.textContent = level.label;
        }
        slider.setAttribute("aria-valuetext", level.label);
        document.querySelectorAll("#mlst-png-res-ticks span").forEach((tick) => {
          tick.classList.toggle("is-active", tick.getAttribute("data-level") === String(idx));
        });
        if (store) {
          try {
            localStorage.setItem(PNG_RES_KEY, String(idx));
          } catch (_err) {}
        }
      };
      slider.addEventListener("input", () => apply(true));
      slider.addEventListener("change", () => apply(true));
      document.querySelectorAll("#mlst-png-res-ticks span").forEach((tick) => {
        tick.addEventListener("click", () => {
          slider.value = tick.getAttribute("data-level") || "1";
          apply(true);
        });
      });
      let saved = "1";
      try {
        saved = localStorage.getItem(PNG_RES_KEY) || "1";
      } catch (_err) {}
      slider.value = saved;
      apply(false);
    }

    const persistDims = () => {
      try {
        localStorage.setItem(
          PNG_DIMS_KEY,
          JSON.stringify({
            width: Number(wEl && wEl.value) || "",
            height: Number(hEl && hEl.value) || "",
          })
        );
      } catch (_err) {}
    };
    [wEl, hEl].forEach((el) => {
      if (!el || el.dataset.bound) {
        return;
      }
      el.dataset.bound = "1";
      el.addEventListener("input", persistDims);
      el.addEventListener("change", persistDims);
    });
    try {
      const stored = JSON.parse(localStorage.getItem(PNG_DIMS_KEY) || "{}") || {};
      if (wEl && stored.width) {
        wEl.value = String(stored.width);
      }
      if (hEl && stored.height) {
        hEl.value = String(stored.height);
      }
    } catch (_err) {}

    const dlSt = document.getElementById("mlst-dl-st");
    const dlAl = document.getElementById("mlst-dl-alleles");
    const dlBoth = document.getElementById("mlst-dl-both");
    if (dlSt && !dlSt.dataset.bound) {
      dlSt.dataset.bound = "1";
      dlSt.addEventListener("click", () => downloadChart("mlst-chart-st", "mlst_st_distribution.png"));
    }
    if (dlAl && !dlAl.dataset.bound) {
      dlAl.dataset.bound = "1";
      dlAl.addEventListener("click", () => downloadChart("mlst-chart-alleles", "mlst_allele_profiles.png"));
    }
    if (dlBoth && !dlBoth.dataset.bound) {
      dlBoth.dataset.bound = "1";
      dlBoth.addEventListener("click", async () => {
        await downloadChart("mlst-chart-st", "mlst_st_distribution.png");
        await downloadChart("mlst-chart-alleles", "mlst_allele_profiles.png");
      });
    }
  }

  function renderCharts(charts) {
    const wrap = document.getElementById("mlst-charts");
    const bar = document.getElementById("mlst-chart-bar");
    const hint = document.getElementById("mlst-charts-hint");
    const stEl = document.getElementById("mlst-chart-st");
    const alEl = document.getElementById("mlst-chart-alleles");
    if (!wrap || !stEl || !alEl || !charts || !window.Plotly) {
      return;
    }
    wrap.style.display = "grid";
    if (hint) {
      hint.style.display = "block";
    }
    if (bar) {
      bar.style.display = "flex";
      const titleEl = bar.querySelector(".bf-plot-bar-title");
      if (titleEl) {
        const n = charts.n_genomes || 0;
        titleEl.textContent =
          "ST distribution (" + n + " genome" + (n === 1 ? "" : "s") + ") and allele profiles";
      }
    }
    bindDownloadControls();

    const colors = plotColors();
    const layoutBase = {
      autosize: true,
      height: CHART_HEIGHT,
      paper_bgcolor: colors.paper,
      plot_bgcolor: colors.paper,
      font: { color: colors.font, size: 11 },
      margin: { t: 40, r: 16, b: 48, l: 48 },
      xaxis: { gridcolor: colors.grid, zerolinecolor: colors.grid },
      yaxis: { gridcolor: colors.grid, zerolinecolor: colors.grid, rangemode: "tozero" },
    };
    const config = { responsive: false, displayModeBar: false, staticPlot: false };

    purge(stEl);
    purge(alEl);
    stEl.dataset.bfNoAutoResize = "1";
    alEl.dataset.bfNoAutoResize = "1";

    if ((charts.st_labels || []).length) {
      const n = charts.n_genomes || 0;
      Plotly.newPlot(
        stEl,
        [{
          type: "bar",
          x: charts.st_labels.map((s) => (s === "-" ? "untyped" : "ST" + s)),
          y: charts.st_counts,
          marker: { color: colors.bar },
        }],
        Object.assign({}, layoutBase, {
          title: "Sequence types (" + n + " genomes)",
        }),
        config
      );
    }

    if ((charts.allele_genes || []).length && (charts.genomes || []).length) {
      Plotly.newPlot(
        alEl,
        [{
          type: "heatmap",
          z: charts.allele_matrix,
          x: charts.allele_genes,
          y: charts.genomes,
          colorscale: "Viridis",
          hoverongaps: false,
        }],
        Object.assign({}, layoutBase, {
          title: "Allele profiles",
          margin: { t: 40, r: 16, b: 80, l: 100 },
        }),
        config
      );
    }
  }

  function render(data) {
    const wrap = document.getElementById("mlstDiv");
    const out = document.getElementById("output-mlst");
    const hint = document.getElementById("mlst-table-hint");
    if (!out || !data || !data.exists) {
      if (wrap) {
        wrap.style.display = "none";
      }
      return;
    }
    if (wrap) {
      wrap.style.display = "block";
    }
    if (hint) {
      const n = data.n_genomes || 0;
      const typed = data.n_typed || 0;
      hint.textContent =
        "PubMLST typing for " +
        n +
        " genome(s); " +
        typed +
        " received a sequence type. Allele columns are scheme-specific.";
    }
    if (data.mlst_table) {
      out.innerHTML = data.mlst_table;
      setTimeout(() => {
        if (!(window.$ && $.fn && $.fn.DataTable)) {
          return;
        }
        if ($.fn.DataTable.isDataTable("#mlst-tab")) {
          $("#mlst-tab").DataTable().destroy();
        }
        $("#mlst-tab").DataTable({
          paging: true,
          pageLength: 10,
          searching: true,
          ordering: true,
          lengthMenu: [
            [10, 25, 50, -1],
            [10, 25, 50, "All"],
          ],
          scrollX: true,
        });
      }, 200);
    } else {
      out.innerHTML = "<p class='text-muted'>MLST finished but no rows were parsed.</p>";
    }
    renderCharts(data.charts);
  }

  window.MlstViz = { render, renderCharts, downloadChart };
})();
