/* MLST table + Plotly ST / allele charts for Post-assembly */
(function () {
  const CHART_HEIGHT = 340;
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
      ? {
          font: "#1f2937",
          grid: "#e5e7eb",
          paper: "#ffffff",
          bar: "#0d9488",
          barLine: "#0f766e",
          muted: "#6b7280",
        }
      : {
          font: "#e5e7eb",
          grid: "#334155",
          paper: "rgba(0,0,0,0)",
          bar: "#2dd4bf",
          barLine: "#14b8a6",
          muted: "#9aa8b5",
        };
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

  function stLabel(st) {
    if (st === "-" || st === "" || st === "0" || st == null) {
      return "Untyped";
    }
    return "ST-" + st;
  }

  function renderStChart(stEl, charts, colors) {
    const labels = (charts.st_labels || []).map(stLabel);
    const counts = charts.st_counts || [];
    if (!labels.length) {
      return;
    }
    const maxCount = Math.max(1, ...counts);
    const few = labels.length <= 2;

    Plotly.newPlot(
      stEl,
      [
        {
          type: "bar",
          x: labels,
          y: counts,
          text: counts.map((c) => String(c) + (c === 1 ? " genome" : " genomes")),
          textposition: "outside",
          cliponaxis: false,
          hovertemplate: "<b>%{x}</b><br>%{y} genome%{customdata}<extra></extra>",
          customdata: counts.map((c) => (c === 1 ? "" : "s")),
          width: few ? Math.min(0.42, 0.85 / Math.max(1, labels.length)) : undefined,
          marker: {
            color: colors.bar,
            line: { color: colors.barLine, width: 1 },
            opacity: 0.92,
          },
        },
      ],
      {
        autosize: true,
        height: CHART_HEIGHT,
        paper_bgcolor: colors.paper,
        plot_bgcolor: colors.paper,
        font: { color: colors.font, size: 11 },
        title: {
          text: "How many genomes share each sequence type",
          font: { size: 13, color: colors.font },
          x: 0.02,
          xanchor: "left",
        },
        margin: { t: 48, r: 20, b: 56, l: 56 },
        bargap: few ? 0.55 : 0.28,
        xaxis: {
          title: { text: "Sequence type (PubMLST)", font: { size: 11, color: colors.muted } },
          type: "category",
          gridcolor: colors.grid,
          zeroline: false,
          tickangle: labels.length > 6 ? -30 : 0,
        },
        yaxis: {
          title: { text: "Genomes", font: { size: 11, color: colors.muted } },
          gridcolor: colors.grid,
          zerolinecolor: colors.grid,
          rangemode: "tozero",
          range: [0, maxCount + Math.max(1, maxCount * 0.35)],
          dtick: maxCount <= 8 ? 1 : undefined,
          tick0: 0,
        },
      },
      { responsive: false, displayModeBar: false, staticPlot: false }
    );
  }

  function showAlleleEmpty(alEl, colors, message) {
    if (!alEl) {
      return;
    }
    alEl.innerHTML =
      "<div class='mlst-allele-empty' style='display:flex;align-items:center;justify-content:center;" +
      "height:100%;min-height:280px;color:" +
      (colors.muted || "#9aa8b5") +
      ";padding:1.25rem;text-align:center;font-size:0.9rem;'>" +
      message +
      "</div>";
  }

  function renderAlleleChart(alEl, charts, colors) {
    const genes = charts.allele_genes || [];
    const genomes = charts.genomes || [];
    let z = charts.allele_matrix || [];
    let text = charts.allele_text || [];
    if (!genomes.length) {
      showAlleleEmpty(alEl, colors, "No MLST genomes to plot.");
      return;
    }
    if (!genes.length) {
      showAlleleEmpty(
        alEl,
        colors,
        "Allele columns were not found in mlst.tsv (only scheme/ST). Re-run MLST so the full allele table is written."
      );
      return;
    }

    // Coerce null/undefined cells so Plotly always has a drawable matrix.
    z = z.map((row, ri) =>
      genes.map((_, ci) => {
        const v = row && row[ci];
        const n = Number(v);
        return Number.isFinite(n) ? n : 0;
      })
    );
    const srcText = text;
    text = z.map((row, ri) =>
      genes.map((_, ci) => {
        const t = srcText[ri] && srcText[ri][ci];
        return t == null || t === "" ? String(row[ci] ?? "-") : String(t);
      })
    );

    const nRows = genomes.length;
    // Pad a single-genome heatmap with blank spacer rows so the real row
    // stays a readable band (Plotly otherwise draws a hairline / empty cell).
    let yLabels = genomes.slice();
    let zPlot = z.map((row) => row.slice());
    let textPlot = text.map((row) => row.slice());
    if (nRows === 1) {
      const blank = genes.map(() => null);
      const blankTxt = genes.map(() => "");
      yLabels = ["", genomes[0], ""];
      zPlot = [blank, z[0], blank];
      textPlot = [blankTxt, text[0], blankTxt];
    }

    Plotly.newPlot(
      alEl,
      [
        {
          type: "heatmap",
          z: zPlot,
          x: genes,
          y: yLabels,
          text: textPlot,
          texttemplate: "%{text}",
          textfont: {
            size: genes.length > 10 ? 9 : 12,
            color: "#f8fafc",
          },
          colorscale: [
            [0, "#0f172a"],
            [0.35, "#134e4a"],
            [0.65, "#14b8a6"],
            [1, "#fde68a"],
          ],
          showscale: true,
          colorbar: {
            title: { text: "Allele #", side: "right" },
            thickness: 12,
            len: 0.75,
            tickfont: { size: 10 },
          },
          hoverongaps: false,
          hovertemplate:
            "<b>%{y}</b><br>Locus: %{x}<br>Allele: %{text}<extra></extra>",
          xgap: 3,
          ygap: 3,
        },
      ],
      {
        autosize: true,
        height: CHART_HEIGHT,
        paper_bgcolor: colors.paper,
        plot_bgcolor: colors.paper,
        font: { color: colors.font, size: 11 },
        title: {
          text: "PubMLST allele IDs at each locus (printed in cells)",
          font: { size: 13, color: colors.font },
          x: 0.02,
          xanchor: "left",
        },
        margin: { t: 48, r: 72, b: 80, l: 120 },
        xaxis: {
          title: { text: "Locus", font: { size: 11, color: colors.muted } },
          type: "category",
          side: "bottom",
          tickangle: genes.length > 8 ? -35 : 0,
          gridcolor: colors.grid,
        },
        yaxis: {
          title: {
            text: nRows === 1 ? "Genome" : "Genomes",
            font: { size: 11, color: colors.muted },
          },
          type: "category",
          categoryorder: "array",
          categoryarray: yLabels,
          autorange: "reversed",
          gridcolor: colors.grid,
          tickfont: { size: 11 },
        },
      },
      { responsive: false, displayModeBar: false, staticPlot: false }
    );
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
      const n = charts.n_genomes || 0;
      const typed = charts.n_typed || 0;
      hint.textContent =
        "Left: count of genomes per sequence type (ST). Right: PubMLST allele numbers at each locus — " +
        n +
        " genome" +
        (n === 1 ? "" : "s") +
        ", " +
        typed +
        " typed.";
    }
    if (bar) {
      bar.style.display = "flex";
      const titleEl = bar.querySelector(".bf-plot-bar-title");
      if (titleEl) {
        const n = charts.n_genomes || 0;
        titleEl.textContent =
          "MLST charts — " + n + " genome" + (n === 1 ? "" : "s");
      }
    }
    bindDownloadControls();

    const colors = plotColors();
    purge(stEl);
    purge(alEl);
    stEl.dataset.bfNoAutoResize = "1";
    alEl.dataset.bfNoAutoResize = "1";

    if ((charts.st_labels || []).length) {
      renderStChart(stEl, charts, colors);
    } else {
      stEl.innerHTML = "";
    }
    // Always attempt allele panel (shows a clear empty-state if no loci).
    renderAlleleChart(alEl, charts, colors);
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
