

const runButton = document.getElementById("run-bt");
const stopButton = document.getElementById("stop-bt");
const helpButton = document.getElementById("help-bt");
let eventSource = null;


const filePicker = document.getElementById("fastq_dir") || document.getElementById("filePicker");
const fastqDirInput = document.getElementById("fastq_dir") || document.getElementById("filePicker");


// check bactflow installation 

function checkBactFlow(){
  fetch('/bactflow-check')
    .then(response => response.json())
    .then(data => {
      const noInstallDiv = document.getElementById("notInstalled");
      const InstalledDiv = document.getElementById("Installed");
      const spinDiv = document.getElementById("installing");
      const instBtn = document.getElementById("install-bt");
      
      console.log(data.status);//for debugging
      if(data.status === "success" && data.installed){
        console.log("✅ BactFlow detected!");//for debugging
        InstalledDiv.innerHTML = `✅ Bactflow is installed`;
        InstalledDiv.style.backgroundColor = "white";
        InstalledDiv.style.display = "flex";
        InstalledDiv.style.color="white";
        InstalledDiv.style.width = "auto";
        InstalledDiv.style.padding = "10px 20px";
        spinDiv.style.display = "none";
        instBtn.style.display = "none";
        noInstallDiv.style.display = "none";
      } else if(data.status === "success" && !data.installed){
        console.log("⏳ BactFlow not installed yet. Checking again...");
        noInstallDiv.style.display = "block";

        setTimeout(checkBactFlow, 2000);

      } 
    })
    .catch(error => {console.log("Error fetching bactflow status:", error);
    setTimeout(checkBactFlow, 2000);
    });
}

document.addEventListener("DOMContentLoaded", function(){
  checkBactFlow();
});

// install button 
const instBtn = document.getElementById('install-bt');
instBtn.addEventListener("click", function(){

  // const outputDiv = document.getElementById("install-output");
  const spinDiv = document.getElementById("installing");
  const form = document.getElementById('preForm');
  const formData = new FormData(form);
  const noInstallDiv = document.getElementById("notInstalled");
  const InstalledDiv = document.getElementById("Installed");
  const instBtn = document.getElementById("install-bt");
  



  fetch("/install-bactflow", {
    method: "POST", 
    body: formData
  })
  .then((response) => response.json())
  .then(data => {
   console.log("data running"+data.running);
   console.log("Server Response:", data); 
    if(data.running){
      spinDiv.style.display = "block";
      instBtn.style.display = "none";
    }
    
   
})
  .catch((error) => {

  outputDiv.innerHTML += "<br>Error: Installation request failed!";
  console.error(error);
})


})

// getting values of the run buttons 
let clickedButtonValue = null;
document.addEventListener("click", function(event) {
  if (event.target.name === "action-assem") {
       clickedButtonValue = event.target.value;
     
  }
});

// trim list
function trimList(){
  const form = document.getElementById('preForm');
  const formData = preFormData();
  const fastqSpin = document.getElementById("spin-fastq");
  const trimBt = document.getElementById("trim-list");
  fastqSpin.style.display = "block";
  trimBt.style.width = "auto";
  trimBt.style.height = "auto";

  fetch("/trim-list", {
    method: "POST",
    body: formData
  })
  .then(response => {
    if(!response.ok){
      return response.json().then(data => {
        throw new Error(data.error);
      });
    }
    return response.json();
  })
  .then(data => {
   
   
    if (data.status == "completed"){
      fastqSpin.style.display = "none";
      // trimBt.style.display = "block";
      trimBt.innerHTML = "Trimming completed";
      trimBt.style.width = "auto";
      trimBt.style.backgroundColor = "green";
      trimBt.style.color = "white";
    } else {  
      trimBt.innerHTML = "Error trimming";
      trimBt.style.width = "auto";
      trimBt.style.backgroundColor = "red";
    }
  }
  )
  .catch(error => {
     console.error("Error trimming list:", error);
    fastqSpin.style.display = "none";
    trimBt.style.display = "block";
    trimBt.innerHTML = "❌ Error Trimming";
    trimBt.classList.remove("btn-warning", "btn-success");
    trimBt.classList.add("btn-danger");
  });
};

function jumpToSection(id) {
  const el = typeof id === "string" ? document.getElementById(id) : id;
  if (!el) {
    return;
  }
  // Results live in a fixed-height overflow column — scroll that pane, not only the window.
  const pane =
    document.getElementById("pre-results-pane") ||
    el.closest(".col-md-9") ||
    null;
  requestAnimationFrame(() => {
    if (pane && pane.scrollHeight > pane.clientHeight) {
      const top =
        el.getBoundingClientRect().top -
        pane.getBoundingClientRect().top +
        pane.scrollTop -
        12;
      pane.scrollTo({ top: Math.max(0, top), behavior: "smooth" });
    } else {
      el.scrollIntoView({ behavior: "smooth", block: "start" });
    }
  });
}

// initiate concatenation
function catFastq(){
  const lsfq = document.getElementById("fastq-ls");
  const lsDiv = document.getElementById("fastq-ls-div");
  const form = document.getElementById('preForm');
  const spinFastqLs = document.getElementById("spin-fastq-list");
  const formData = preFormData();
  spinFastqLs.style.display = "flex";
  jumpToSection("spin-fastq-list");
  fetch("/ls-fastq", {
    method: "POST",
    body: formData
  })
  .then(response => {
    if(!response.ok){
      return response.json().then(data => {
        throw new Error(data.error);
      });
    }
    return response.json();
  })
  .then(data => {
    
    if (data.html_table){
      spinFastqLs.style.display = "none";
      lsDiv.style.display = "block";
      lsfq.innerHTML = data.html_table;
      if (data.note) {
        const note = document.createElement("p");
        note.style.marginTop = "8px";
        note.innerHTML = data.note;
        lsfq.prepend(note);
      }
      if (data.looks_like_subreads) {
        const kind = document.getElementById("pacbio_read_kind");
        const platform = document.getElementById("read_platform");
        if (kind) {
          kind.value = "clr";
        }
        if (platform && platform.value === "auto") {
          platform.value = "pacbio";
          updatePlatformUI();
        }
      }
      if (data.extension) {
        const extInput = document.getElementById("extension");
        if (extInput && (extInput.value === ".fastq.gz" || !extInput.value)) {
          extInput.value = data.extension;
        }
      }
      jumpToSection("fastq-ls-div");
      setTimeout(() => {
        if ($.fn.DataTable) {
          console.log("✅ DataTables is loaded, initializing table...");
          $('#read-fastq').DataTable({
            "paging": true,
            "pageLength": 10,
            "searching": true,
            "ordering": true,
            "lengthMenu": [[10, 25, 50, -1], [10, 25, 50, "All"]],
            "responsive": true
          }); } else {
            console.error("❌ DataTables is not loaded, retrying in 2ms...");
        }
      }, 20); 
    } else {
      lsfq.innerHTML = `<p>Error: ${data.error}</p>`;
      lsDiv.style.display = "block";
      jumpToSection("fastq-ls-div");
    }

    document.getElementById("size-slider").addEventListener("input", function() {
      let threshold = parseInt(this.value);
      document.getElementById("size-value").textContent = threshold;
      
      document.querySelectorAll("#read-fastq tbody tr").forEach((row, index) => {
   
          // if (index === 0) return; // Skip header row
          let sizeFq = row.cells[2].textContent.trim();
          let sizeMb = convertToMB(sizeFq);
          console.log(`Size: ${sizeFq} to ${sizeMb} MB`);
          if (sizeMb < threshold) {
              row.classList.add("low-seqs");
          } else {
              row.classList.remove("low-seqs");
          }
      });
  });
  
   
  })
  .catch(error => {
    console.error("Error fetching fastq files:", error);
    spinFastqLs.style.display = "none";
    lsDiv.style.display = "block";
    lsfq.innerHTML = `<p style="color:#b00020;"><strong>Could not list FASTQ files.</strong><br>${error.message}</p>`;
    jumpToSection("fastq-ls-div");
  });
}

function convertToMB(sizeText) { 
  let match = sizeText.match(/([\d.]+)\s*([MGK]?)/i); // Extract number & unit
  if (!match) return 0;

  let size = parseFloat(match[1]); // Extract numeric part
  let unit = match[2].toUpperCase(); // Extract unit (M, G, K)

  if (unit === "G") size *= 1024; // Convert GB to MB
  if (unit === "K") size /= 1024;  // Convert KB to MB

  return size; // Return MB
}

function preFormData() {
  const form = document.getElementById("preForm");
  const formData = new FormData(form);
  if ((formData.get("read_platform") || "") === "illumina") {
    formData.set("concat_reads", "false");
  }
  return formData;
}

let lastReadPlatform = null;

function destroyTable(selector) {
  if (window.$ && $.fn.DataTable && $(selector).length && $.fn.DataTable.isDataTable(selector)) {
    $(selector).DataTable().destroy();
  }
}

function purgePlot(el) {
  if (!el) {
    return;
  }
  if (typeof Plotly !== "undefined" && typeof Plotly.purge === "function") {
    try {
      Plotly.purge(el);
    } catch (err) {
      // ignore stale plot handles
    }
  }
  el.style.display = "none";
}

function clearPlots() {
  ["qual-plot", "qual-heat-plot", "qual-heat-pool-plot"].forEach((id) => {
    purgePlot(document.getElementById(id));
  });
  ["bar-qual-plot", "bar-qual-heat-plot", "bar-qual-heat-pool-plot",
    "range1Div", "range2Div", "range3Div"].forEach((id) => {
    const el = document.getElementById(id);
    if (el) {
      el.style.display = "none";
    }
  });
}

function clearResultsPane() {
  ["spin-vis", "spin-seqkit", "spin-fastq-list", "spin-fastq"].forEach((id) => {
    const el = document.getElementById(id);
    if (el) {
      el.style.display = "none";
    }
  });
  const visError = document.getElementById("vis-error");
  if (visError) {
    visError.style.display = "none";
    visError.textContent = "";
    visError.innerHTML = "";
  }
  destroyTable("#stats-table");
  destroyTable("#read-fastq");
  const seqDiv = document.getElementById("seq-div");
  const outputDiv = document.getElementById("output-seqkit");
  const threshIn = document.getElementById("threhold-input");
  if (seqDiv) {
    seqDiv.style.display = "none";
  }
  if (outputDiv) {
    outputDiv.innerHTML = "";
  }
  if (threshIn) {
    threshIn.style.display = "none";
  }
  const lsDiv = document.getElementById("fastq-ls-div");
  const lsfq = document.getElementById("fastq-ls");
  if (lsDiv) {
    lsDiv.style.display = "none";
  }
  if (lsfq) {
    lsfq.innerHTML = "<h4 class=\"text-left mb-3\">List of fastq files <strong>(FASTQ)</strong></h4>";
  }
  clearPlots();
  const trimBt = document.getElementById("trim-list");
  if (trimBt) {
    trimBt.innerHTML = "<strong>Trim the list!</strong>";
    trimBt.style.backgroundColor = "";
    trimBt.style.color = "";
    trimBt.classList.remove("btn-danger", "btn-success");
    trimBt.classList.add("btn-warning");
  }
}

function selectedPlatform() {
  return ((document.getElementById("read_platform") || {}).value || "auto").toLowerCase();
}

function selectedPacbioKind() {
  return ((document.getElementById("pacbio_read_kind") || {}).value || "hifi").toLowerCase();
}

function resolvePlotPlatform(declared) {
  const p = (declared || selectedPlatform() || "auto").toLowerCase();
  if (p === "illumina") {
    return "illumina";
  }
  if (p === "pacbio") {
    return "pacbio";
  }
  return "ont";
}

function updatePlatformUI() {
  const platform = selectedPlatform();
  const illumina = platform === "illumina";
  const pacbio = platform === "pacbio";
  const concatGroup = document.getElementById("concatReadsGroup");
  const concatSel = document.getElementById("concat_reads");
  const catBtn = document.getElementById("read-cat");
  const filterBtn = document.getElementById("illumina-filter-bt");
  const filterGroup = document.getElementById("illuminaFilterGroup");
  const filterParams = document.getElementById("illuminaFilterParams");
  const pacbioKindGroup = document.getElementById("pacbioKindGroup");
  const filterOn = (document.getElementById("illumina_filter") || {}).value === "true";

  if (concatGroup) {
    concatGroup.style.display = illumina ? "none" : "";
  }
  if (concatSel && illumina) {
    concatSel.value = "false";
  }
  if (concatSel && lastReadPlatform === "illumina" && !illumina) {
    concatSel.value = "true";
  }
  if (concatSel && pacbio) {
    concatSel.value = "true";
  }
  if (catBtn) {
    catBtn.style.display = illumina ? "none" : "";
  }
  if (filterBtn) {
    filterBtn.style.display = illumina ? "" : "none";
  }
  if (filterGroup) {
    filterGroup.style.display = illumina ? "" : "none";
  }
  if (filterParams) {
    filterParams.style.display = illumina && filterOn ? "" : "none";
  }
  if (pacbioKindGroup) {
    pacbioKindGroup.style.display = pacbio ? "" : "none";
  }
  if (lastReadPlatform !== null && lastReadPlatform !== platform) {
    clearResultsPane();
    setPlotChrome(resolvePlotPlatform(platform));
  }
  lastReadPlatform = platform;
}

function showNotice(msg, isError) {
  const el = document.getElementById("vis-error");
  if (!el) {
    return;
  }
  el.style.display = "block";
  el.style.background = isError ? "#fdeeee" : "#e8f7f4";
  el.style.color = isError ? "#8b1e1e" : "#0e5f58";
  el.style.borderColor = isError ? "#e0a0a0" : "#12b5a8";
  el.innerHTML = msg;
}

function filterIllumina() {
  const btn = document.getElementById("illumina-filter-bt");
  const filterSel = document.getElementById("illumina_filter");
  if (filterSel) {
    filterSel.value = "true";
  }
  updatePlatformUI();
  if (btn) {
    btn.disabled = true;
    btn.textContent = "Filtering…";
  }
  fetch("/filter-illumina", {
    method: "POST",
    body: preFormData()
  })
    .then(parseJsonResponse)
    .then((data) => {
      if (btn) {
        btn.disabled = false;
        btn.innerHTML = "<strong>Filter Illumina reads</strong>";
      }
      showNotice(data.message || "Illumina reads filtered.", false);
    })
    .catch((error) => {
      if (btn) {
        btn.disabled = false;
        btn.innerHTML = "<strong>Filter Illumina reads</strong>";
      }
      showNotice(error.message || String(error), true);
    });
}

function parseJsonResponse(response) {
  return response.text().then((text) => {
    let data = {};
    if (text) {
      try {
        data = JSON.parse(text);
      } catch (err) {
        throw new Error(text.replace(/<[^>]+>/g, " ").slice(0, 400) || `HTTP ${response.status}`);
      }
    }
    if (!response.ok) {
      throw new Error(data.error || data.message || `HTTP ${response.status}`);
    }
    return data;
  });
}

function showSeqError(msg) {
  const spinSeq = document.getElementById("spin-seqkit");
  const seqDiv = document.getElementById("seq-div");
  const outputDiv = document.getElementById("output-seqkit");
  if (spinSeq) {
    spinSeq.style.display = "none";
  }
  if (seqDiv) {
    seqDiv.style.display = "block";
  }
  if (outputDiv) {
    outputDiv.innerHTML = `<p style="color:#b00020;"><strong>Could not generate read statistics.</strong><br>${msg}</p>`;
  }
}

function downloadPlot(plotId, filename, format) {
  const el = document.getElementById(plotId);
  if (!el || typeof Plotly === "undefined" || typeof Plotly.downloadImage !== "function") {
    return;
  }
  const kind = format === "svg" ? "svg" : "png";
  Plotly.downloadImage(el, {
    format: kind,
    filename: filename,
    height: Math.max(el.clientHeight || 800, 600),
    width: Math.max(el.clientWidth || 1200, 800),
    scale: kind === "png" ? 2 : 1
  });
}

function setPlotlyTitle(plotData, text) {
  if (!plotData) {
    return plotData;
  }
  if (!plotData.layout) {
    plotData.layout = {};
  }
  const prev = plotData.layout.title;
  if (prev && typeof prev === "object") {
    plotData.layout.title = Object.assign({}, prev, { text: text });
  } else {
    plotData.layout.title = { text: text };
  }
  if (Array.isArray(plotData.layout.annotations)) {
    plotData.layout.annotations = plotData.layout.annotations.map((ann) => {
      if (ann && typeof ann.text === "string" && /\bONT\b/.test(ann.text)) {
        return Object.assign({}, ann, { text: ann.text.replace(/\bONT\b/g, "PacBio") });
      }
      return ann;
    });
  }
  return plotData;
}

function visualizationPlatform(dataPlatform) {
  if (selectedPlatform() === "pacbio" || (dataPlatform || "").toLowerCase() === "pacbio") {
    return "pacbio";
  }
  return resolvePlotPlatform(dataPlatform || selectedPlatform());
}

function setPlotChrome(platform, pacbioKind, pacbioVizMode) {
  const p = resolvePlotPlatform(platform);
  let specs;
  if (p === "illumina") {
    specs = [
      ["bar-qual-plot", "qual-plot", "Illumina per-cycle quality", "illumina_per_cycle_quality"],
      ["bar-qual-heat-plot", "qual-heat-plot", "Illumina mean quality distribution", "illumina_quality_histogram"],
      ["bar-qual-heat-pool-plot", "qual-heat-pool-plot", "Illumina per-cycle nucleotide content", "illumina_base_composition"]
    ];
  } else if (p === "pacbio") {
    specs = [
      ["bar-qual-plot", "qual-plot", "PacBio read quality box plot", "pacbio_read_quality_boxplot"],
      ["bar-qual-heat-plot", "qual-heat-plot", "PacBio read length vs quality (per sample)", "pacbio_length_quality_heatmap"],
      ["bar-qual-heat-pool-plot", "qual-heat-pool-plot", "PacBio read length vs quality (pooled)", "pacbio_length_quality_pooled"]
    ];
  } else {
    specs = [
      ["bar-qual-plot", "qual-plot", "ONT read quality box plot", "ont_read_quality_boxplot"],
      ["bar-qual-heat-plot", "qual-heat-plot", "ONT read length vs quality (per sample)", "ont_length_quality_heatmap"],
      ["bar-qual-heat-pool-plot", "qual-heat-pool-plot", "ONT read length vs quality (pooled)", "ont_length_quality_pooled"]
    ];
  }
  specs.forEach(([barId, plotId, title, fname]) => {
    const bar = document.getElementById(barId);
    if (!bar) {
      return;
    }
    const titleEl = bar.querySelector(".bf-plot-bar-title");
    if (titleEl) {
      titleEl.textContent = title;
    }
    const buttons = bar.querySelectorAll("button");
    if (buttons[0]) {
      buttons[0].setAttribute("onclick", `downloadPlot('${plotId}','${fname}','png')`);
    }
    if (buttons[1]) {
      buttons[1].setAttribute("onclick", `downloadPlot('${plotId}','${fname}','svg')`);
    }
  });
}

function read_stat(){
  const spinSeq = document.getElementById("spin-seqkit");
  spinSeq.style.display = "flex";
  jumpToSection("spin-seqkit");
  const form = document.getElementById('preForm');
  const formData = preFormData();
  fetch("/reads-stat", {
    method: "POST",
    body: formData
  })
    .then(parseJsonResponse)
    .then(data => {
      const seqDiv = document.getElementById("seq-div");
      const outputDiv = document.getElementById("output-seqkit");
      const threshIn = document.getElementById("threhold-input");
      if (!data.html_output) {
        throw new Error(data.error || data.message || "Stats generation failed");
      }
      spinSeq.style.display = "none";
      seqDiv.style.display = "block";
      threshIn.style.display = "block";
      outputDiv.innerHTML = data.html_output;
      jumpToSection("seq-div");

      let columnIndex = -1;
      document.querySelectorAll("#stats-table thead th").forEach((th, index) => {
          if (th.textContent.trim() === "sum_len") {
              columnIndex = index;
          }
      });
       
      if (columnIndex !== -1) {
        document.getElementById("depth-slider").addEventListener("input", function() {
          let threshold = parseInt(this.value);
          document.getElementById("depth-value").textContent = threshold;
          document.querySelectorAll("#stats-table tbody tr").forEach(row => {
              let numSeqs = parseInt(row.cells[columnIndex].textContent.replace(/,/g, '')) || 0;
              if (numSeqs < threshold) {
                  row.classList.add("low-seqs");
              } else {
                  row.classList.remove("low-seqs");
              }
          });
        });
      }

      setTimeout(() => {
        if ($.fn.DataTable) {
          if ($.fn.DataTable.isDataTable("#stats-table")) {
            $("#stats-table").DataTable().destroy();
          }
          $('#stats-table').DataTable({
            "paging": true,
            "pageLength": 10,
            "searching": true,
            "ordering": true,
            "lengthMenu": [[10, 25, 50, -1], [10, 25, 50, "All"]],
            "responsive": true
          });
        }
        jumpToSection("seq-div");
      }, 20);
    })
    .catch(error => {
      console.error("Error fetching read-stat:", error);
      showSeqError(error.message || String(error));
      jumpToSection("seq-div");
    });
}


// read visualization 

function read_vis(){
  const qualDiv = document.getElementById("qual-plot");
  const form = document.getElementById("preForm");
  const formData = preFormData();
  const spinVis = document.getElementById("spin-vis");
  const visError = document.getElementById("vis-error");
  const qualHeatDiv = document.getElementById("qual-heat-plot");
  const qualHeatDivPool = document.getElementById("qual-heat-pool-plot");
  const range1Div = document.getElementById("range1Div");
  const range2Div = document.getElementById("range2Div");
  const range3Div = document.getElementById("range3Div");
  const plotBars = [
    document.getElementById("bar-qual-plot"),
    document.getElementById("bar-qual-heat-plot"),
    document.getElementById("bar-qual-heat-pool-plot")
  ];

  clearPlots();
    setPlotChrome(visualizationPlatform(selectedPlatform()), selectedPacbioKind());
  if (visError) {
    visError.style.display = "none";
    visError.textContent = "";
    visError.innerHTML = "";
  }
  spinVis.style.display = "flex";
  jumpToSection("spin-vis");

  fetch("/plot-qual", {method: "POST", body: formData})
  .then(parseJsonResponse)
  .then( data =>{
    if (data.error){
      throw new Error(data.error);
    }
    if(!data.graph){
      throw new Error("No plot data returned");
    }
    setPlotChrome(
      visualizationPlatform(data.platform),
      data.pacbio_read_kind || selectedPacbioKind(),
      data.pacbio_viz_mode
    );
    if (data.note) {
      const note = String(data.note);
      const isSub = /subread|CLR|not HiFi/i.test(note);
      showNotice(note, isSub);
    }
    spinVis.style.display = "none";
    qualDiv.style.display = "block";
    qualHeatDiv.style.display = "block";
    qualHeatDivPool.style.display = "block";
    range1Div.style.display = "flex";
    range2Div.style.display = "flex";
    range3Div.style.display = "flex";
    plotBars.forEach((bar) => {
      if (bar) {
        bar.style.display = "flex";
      }
    });

    let width1 = qualDiv.clientWidth;
    let width2 = qualHeatDiv.clientWidth;
    let width3 = qualHeatDivPool.clientWidth;
    let height1 = qualDiv.clientHeight;
    let height2 = qualHeatDiv.clientHeight;
    let height3 = qualHeatDivPool.clientHeight;

    document.getElementById("range1").oninput = function() {
      width1 = parseInt(this.value);
      document.getElementById("width1-value").textContent = width1;
      Plotly.relayout(qualDiv, {width: width1});
    };
    document.getElementById("range1H").oninput = function() {
      height1 = parseInt(this.value);
      document.getElementById("height1-value").textContent = height1;
      Plotly.relayout(qualDiv, {height: height1});
    };
    document.getElementById("range2").oninput = function() {
      width2 = parseInt(this.value);
      document.getElementById("width2-value").textContent = width2;
      Plotly.relayout(qualHeatDiv, {width: width2});
    };
    document.getElementById("range2H").oninput = function() {
      height2 = parseInt(this.value);
      document.getElementById("height2-value").textContent = height2;
      Plotly.relayout(qualHeatDiv, { height: height2});
    };
    document.getElementById("range3").oninput = function() {
      width3 = parseInt(this.value);
      document.getElementById("width3-value").textContent = width3;
      Plotly.relayout(qualHeatDivPool, {width: width3});
    };
    document.getElementById("range3H").oninput = function() {
      height3 = parseInt(this.value);
      document.getElementById("height3-value").textContent = height3;
      Plotly.relayout(qualHeatDivPool, { height: height3});
    };

    var layout1 = { width: width1, height: height1 };
    var layout2 = { width: width1, height: height2 };
    var layout3 = { width: width3, height: height3 };

    var plotData1 = JSON.parse(data.graph);
    var plotData2 = JSON.parse(data.graph_heat);
    var plotData3 = JSON.parse(data.graph_heat2);

    plotData1.layout = { ...plotData1.layout, ...layout1 };
    plotData2.layout = { ...plotData2.layout, ...layout2 };
    plotData3.layout = { ...plotData3.layout, ...layout3 };

    if (visualizationPlatform(data.platform) === "pacbio") {
      setPlotlyTitle(plotData1, "PacBio read quality box plot");
      setPlotlyTitle(plotData2, "PacBio read length vs quality | per sample");
      setPlotlyTitle(plotData3, "PacBio read length vs quality | pooled");
    }

    Plotly.newPlot(qualDiv, plotData1.data, plotData1.layout);
    Plotly.newPlot(qualHeatDiv, plotData2.data, plotData2.layout);
    Plotly.newPlot(qualHeatDivPool, plotData3.data, plotData3.layout);
    const noteVisible = visError && visError.style.display !== "none" && visError.innerHTML;
    jumpToSection(noteVisible ? "vis-error" : "bar-qual-plot");
  })
  .catch(error => {
    console.error("Error fetching plot", error);
    spinVis.style.display = "none";
    if (visError) {
      visError.style.display = "block";
      visError.innerHTML = `<strong>Could not create plots.</strong><br>${error.message || error}`;
      jumpToSection("vis-error");
    }
  });
}

// function update buttons
function updateButtonStates(status) {
  if (status === "running") {
    runButton.disabled = true;
    stopButton.disabled = false;
    helpButton.disabled = true; // Optional: disable help while running
  } else if (status === "stopped") {
    runButton.disabled = false;
    stopButton.disabled = true;
    helpButton.disabled = false;
  } 
}

document.getElementById("preForm").addEventListener("submit", (e) => {
 e.preventDefault();
 outputDiv.innerHTML = "";

  const action = e.submitter.value;

  document.getElementById('run-bt').disabled = true;
  document.getElementById('stop-bt').disabled = false; 
  document.getElementById('help-bt').disabled = true;

  run_wf(action);
});

// saving inputs in the local storage to prevent refereshing
function restoreFormData(){
  let savedData = localStorage.getItem("bactflowFormData");
  if(savedData){
    let formData = JSON.parse(savedData);
    let form = document.getElementById("preForm");
    
    Array.from(form.elemnts).forEach(function (element){
      if (element.name && formData[element.name] !== undefined ) {
        if (elemnt.type === "checkbox" || elemnt.type === "radio"){
          elemnt.checked = formData[element.name];
        } else {
          element.value = formData[element.name];
        }
      }
    });
  }
}



function saveFormData(){
  let form = document.getElementById("preForm");
  let formData = {};

  Array.from(form.elemnts).forEach(function (element){
    if (element.name) {
      if(element.name === "checkbox" || element.type === "radio"){
        formData[element.name] = element.checked;
      } else {
        formData[element.name] = element.value;
      }
    }
  });
  localStorage.setItem("bactflowFormData", JSON.stringify(formData));
}

function clearFormData(){
  localStorage.removeItem("bactflowFormData");
}


document.addEventListener("DOMContentLoaded", function(){
  restoreFormData(); //from local storage
  updatePlatformUI();
  const platform = document.getElementById("read_platform");
  const illuminaFilter = document.getElementById("illumina_filter");
  const pacbioKind = document.getElementById("pacbio_read_kind");
  if (platform) {
    platform.addEventListener("change", updatePlatformUI);
  }
  if (illuminaFilter) {
    illuminaFilter.addEventListener("change", updatePlatformUI);
  }
  if (pacbioKind) {
    pacbioKind.addEventListener("change", () => {
      if (selectedPlatform() === "pacbio") {
        setPlotChrome("pacbio");
      }
    });
  }

  document.getElementById("preForm").addEventListener("input", () =>{
    saveFormData();
  });
});


// Handle directory selection
if (filePicker) {
  filePicker.addEventListener("change", (event) => {
    const selectedFiles = Array.from(event.target.files || []);
    if (selectedFiles.length > 0 && fastqDirInput) {
      const absolutePath = selectedFiles[0].path || selectedFiles[0].webkitRelativePath.split("/")[0];
      fastqDirInput.value = absolutePath;
    }
  });
}


// tggler function
const toggler = (motherId, childId) => {
  document.addEventListener("DOMContentLoaded", ()=>{
  const motherDiv = document.getElementById(motherId);
  const childDiv  = document.getElementById(childId);

  const  toggleDiv = () => {
    if(motherDiv.value === "true"){
      childDiv.style.display = "block";
      childDiv.removeAttribute("disabled"); //so that it still remebers user input
    } else {
      childDiv.style.display = "none";
      childDiv.setAttribute("disabled", "true");
    }
  }

  toggleDiv();

  motherDiv.addEventListener("change", toggleDiv)
});
};


