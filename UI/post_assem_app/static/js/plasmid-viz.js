/**
 * Shared geNomad plasmid + virus visualization (tables + Plotly charts).
 * Used by assembly and post-assembly UIs.
 */
(function (global) {
  function el(id) {
    return document.getElementById(id);
  }

  function hasPlotly() {
    return typeof global.Plotly !== "undefined";
  }

  function initDataTable(tableId) {
    setTimeout(() => {
      if (!(global.$ && $.fn && $.fn.DataTable)) {
        return;
      }
      const sel = "#" + tableId;
      if ($.fn.DataTable.isDataTable(sel)) {
        $(sel).DataTable().destroy();
      }
      $(sel).DataTable({
        paging: true,
        pageLength: 10,
        searching: true,
        ordering: true,
        lengthMenu: [
          [10, 25, 50, -1],
          [10, 25, 50, "All"],
        ],
      });
    }, 300);
  }

  function renderSummary(boxId, summary, kind) {
    const box = el(boxId);
    if (!box || !summary) {
      return;
    }
    const cards =
      kind === "virus"
        ? [
            ["Virus contigs", summary.n_viruses ?? summary.n_items ?? 0],
            ["Genomes with viruses", summary.n_genomes ?? 0],
            ["Mean score", summary.mean_score == null ? "—" : summary.mean_score],
            ["Total virus bp", summary.total_bp ?? 0],
            ["With taxonomy", summary.with_taxonomy ?? 0],
          ]
        : [
            ["Plasmid contigs", summary.n_plasmids ?? summary.n_items ?? 0],
            ["Genomes with plasmids", summary.n_genomes ?? 0],
            ["Mean score", summary.mean_score == null ? "—" : summary.mean_score],
            ["Total plasmid bp", summary.total_bp ?? 0],
            ["With AMR genes", summary.with_amr ?? 0],
            ["With conjugation genes", summary.with_conjugation ?? 0],
          ];
    box.innerHTML = cards
      .map(
        ([label, value]) =>
          `<div class="plasmid-stat"><div class="plasmid-stat-value">${value}</div><div class="plasmid-stat-label">${label}</div></div>`
      )
      .join("");
  }

  function renderCharts(prefix, charts, kind) {
    if (!hasPlotly() || !charts) {
      return;
    }
    const layoutBase = {
      margin: { t: 40, r: 20, b: 60, l: 60 },
      paper_bgcolor: "#0f1419",
      plot_bgcolor: "#0f1419",
      font: { size: 12, color: "#e7edf4", family: "IBM Plex Sans, sans-serif" },
      xaxis: { gridcolor: "#2a3441", zerolinecolor: "#2a3441", color: "#9aa8b5" },
      yaxis: { gridcolor: "#2a3441", zerolinecolor: "#2a3441", color: "#9aa8b5" },
    };
    function layout(extra) {
      const e = extra || {};
      return Object.assign({}, layoutBase, e, {
        xaxis: Object.assign({}, layoutBase.xaxis, e.xaxis || {}),
        yaxis: Object.assign({}, layoutBase.yaxis, e.yaxis || {}),
      });
    }
    const config = { responsive: true, displayModeBar: false };
    const label = kind === "virus" ? "Virus" : "Plasmid";
    const scoreTitle =
      kind === "virus" ? "geNomad virus score distribution" : "geNomad plasmid score distribution";
    const scoreAxis = kind === "virus" ? "Virus score (0–1)" : "Plasmid score (0–1)";
    const colors =
      kind === "virus"
        ? { bar: "#7b2cbf", hist: "#ff6b6b", len: "#4c6ef5", pie: ["#5f3dc4", "#7950f2", "#9775fa", "#b197fc", "#d0bfff"] }
        : { bar: "#2a9d8f", hist: "#e76f51", len: "#457b9d", pie: ["#264653", "#2a9d8f", "#e9c46a", "#f4a261", "#e76f51"] };

    const byGenome = el(prefix + "-chart-by-genome");
    if (byGenome && charts.by_genome && charts.by_genome.genomes.length) {
      Plotly.newPlot(
        byGenome,
        [
          {
            type: "bar",
            x: charts.by_genome.genomes,
            y: charts.by_genome.counts,
            marker: { color: colors.bar },
            name: label + " contigs",
          },
        ],
        layout({
          title: label + " contigs per genome",
          yaxis: { title: "Count", rangemode: "tozero" },
          xaxis: { tickangle: -30 },
        }),
        config
      );
    }

    const scores = el(prefix + "-chart-scores");
    if (scores && charts.scores && charts.scores.values.length) {
      Plotly.newPlot(
        scores,
        [
          {
            type: "histogram",
            x: charts.scores.values,
            nbinsx: Math.min(20, Math.max(5, charts.scores.values.length)),
            marker: { color: colors.hist },
            name: "Score",
          },
        ],
        layout({
          title: scoreTitle,
          xaxis: { title: scoreAxis, range: [0, 1.05] },
          yaxis: { title: "Contigs", rangemode: "tozero" },
        }),
        config
      );
    }

    const lengths = el(prefix + "-chart-lengths");
    if (lengths && charts.lengths && charts.lengths.contigs.length) {
      Plotly.newPlot(
        lengths,
        [
          {
            type: "bar",
            orientation: "h",
            y: charts.lengths.contigs.slice().reverse(),
            x: charts.lengths.lengths.slice().reverse(),
            text: charts.lengths.genomes.slice().reverse(),
            marker: { color: colors.len },
            hovertemplate: "%{y}<br>%{x} bp<br>%{text}<extra></extra>",
          },
        ],
        layout({
          title: "Longest predicted " + label.toLowerCase() + "s (top 40)",
          xaxis: { title: "Length (bp)" },
          yaxis: { automargin: true },
          height: Math.max(320, charts.lengths.contigs.length * 18),
          margin: { t: 40, r: 20, b: 50, l: 120 },
        }),
        config
      );
    }

    const topologyOrTax = el(prefix + "-chart-topology");
    if (topologyOrTax) {
      const useTax = kind === "virus" && charts.taxonomy && charts.taxonomy.labels.length;
      const pie = useTax ? charts.taxonomy : charts.topology;
      const pieTitle = useTax ? "Virus taxonomy" : "Topology";
      if (pie && pie.labels && pie.labels.length) {
        Plotly.newPlot(
          topologyOrTax,
          [
            {
              type: "pie",
              labels: pie.labels,
              values: pie.counts,
              hole: 0.35,
              marker: { colors: colors.pie },
            },
          ],
          layout({
            title: pieTitle,
            margin: { t: 40, r: 10, b: 10, l: 10 },
            showlegend: true,
          }),
          config
        );
      }
    }
  }

  function emptyMessage(kind, genomes) {
    const list = Array.isArray(genomes) ? genomes.filter(Boolean) : [];
    const genomeBit = list.length
      ? list.map((g) => `<code>${g}</code>`).join(", ")
      : "the provided genome(s)";
    if (kind === "virus") {
      return `No viruses / proviruses detected for ${genomeBit}.`;
    }
    return `No plasmids detected for ${genomeBit}.`;
  }

  function renderSection(opts) {
    const {
      sectionId,
      hintId,
      summaryId,
      chartPrefix,
      outId,
      tableHtml,
      charts,
      count,
      kind,
      fastaHint,
      empty,
      genomes,
    } = opts;
    const section = el(sectionId);
    const out = el(outId);
    const hint = el(hintId);
    if (!out) {
      return false;
    }
    if (!tableHtml && !empty) {
      if (section) {
        section.style.display = "none";
      }
      return false;
    }
    if (section) {
      section.style.display = "block";
    }
    const summaryBox = el(summaryId);
    const firstChart = el(chartPrefix + "-chart-by-genome");
    const chartsBox = firstChart ? firstChart.parentElement : null;
    const setShown = (shown) => {
      [summaryBox, chartsBox, out].forEach((node) => {
        if (node) {
          node.style.display = shown ? "" : "none";
        }
      });
    };
    if (empty) {
      setShown(false);
      if (summaryBox) {
        summaryBox.innerHTML = "";
      }
      [
        chartPrefix + "-chart-by-genome",
        chartPrefix + "-chart-scores",
        chartPrefix + "-chart-lengths",
        chartPrefix + "-chart-topology",
      ].forEach((id) => {
        const node = el(id);
        if (node) {
          node.innerHTML = "";
        }
      });
      out.innerHTML = "";
      if (hint) {
        hint.innerHTML = emptyMessage(kind, genomes);
        hint.classList.add("bf-empty-note");
      }
      return true;
    }
    setShown(true);
    if (hint) {
      hint.classList.remove("bf-empty-note");
    }
    if (hint) {
      hint.textContent = count
        ? `geNomad ${kind} calls for ${count} contig(s). ${fastaHint}`
        : `geNomad ${kind} calls.`;
    }
    renderSummary(summaryId, charts && charts.summary, kind);
    renderCharts(chartPrefix, charts, kind);
    out.innerHTML = tableHtml;
    initDataTable(kind === "virus" ? "virus-tab" : "plasmid-tab");
    return true;
  }

  function render(data) {
    const wrap = el("plasmidDiv");
    if (!data || !data.exists) {
      if (wrap) {
        wrap.style.display = "none";
      }
      return;
    }
    if (wrap) {
      wrap.style.display = "block";
    }

    const genomes = data.genomes || [];
    const hasPlasmid = renderSection({
      sectionId: "plasmid-section",
      hintId: "plasmid-table-hint",
      summaryId: "plasmid-summary",
      chartPrefix: "plasmid",
      outId: "output-plasmids",
      tableHtml: data.plasmid_table,
      charts: data.charts,
      count: data.n_plasmids,
      kind: "plasmid",
      fastaHint: "FASTA sequences are in plasmid_out/plasmids.",
      empty: Boolean(data.plasmid_empty),
      genomes,
    });

    const hasVirus = renderSection({
      sectionId: "virus-section",
      hintId: "virus-table-hint",
      summaryId: "virus-summary",
      chartPrefix: "virus",
      outId: "output-viruses",
      tableHtml: data.virus_table,
      charts: data.virus_charts,
      count: data.n_viruses,
      kind: "virus",
      fastaHint: "FASTA sequences are in plasmid_out/viruses.",
      empty: Boolean(data.virus_empty),
      genomes,
    });

    if (!hasPlasmid && !hasVirus && wrap) {
      wrap.style.display = "none";
    }
  }

  global.PlasmidViz = { render, renderCharts, renderSummary };
})(window);
