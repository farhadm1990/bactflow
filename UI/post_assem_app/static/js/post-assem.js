const outputDiv = document.getElementById("output-bactflow");
const runButton = document.getElementById("run-bt");
const stopButton = document.getElementById("stop-bt");
const helpButton = document.getElementById("help-bt");
let eventSource = null;
const STOP_USER_MSG = "Bactflow has been stopped by the user!";
// const form = document.getElementById("postForm");
// const formData = new FormData(form);
const filePicker = document.getElementById("filePicker");
const fastqDirInput = document.getElementById("fastq_dir");






// connect to stream when on assembly
function connectToStream(action){
  // Always (re)connect — a stale EventSource left old process lines on screen.
  if (eventSource){
    try {
      eventSource.onmessage = null;
      eventSource.onerror = null;
      eventSource.close();
    } catch (e) { /* ignore */ }
    eventSource = null;
  }

  BactflowTerminal.init();
  let logBuffer = [];
  let flushTimer = null;
  if (!window.__bactflowStreamGen) {
    window.__bactflowStreamGen = 0;
  }
  const gen = ++window.__bactflowStreamGen;

  const flushLog = () => {
    flushTimer = null;
    if (gen !== window.__bactflowStreamGen) {
      logBuffer = [];
      return;
    }
    if (!logBuffer.length) {
      return;
    }
    BactflowTerminal.appendMany(logBuffer);
    logBuffer = [];
  };

  eventSource = new EventSource(`/stream_bactflow?action-assem=${action}`);
  if (window.BactflowMeters) {
    BactflowMeters.attachStream(eventSource);
  }
  
  eventSource.onmessage = (event) =>{
    if (gen !== window.__bactflowStreamGen) {
      return;
    }
    logBuffer.push(event.data);
    if (!flushTimer) {
      flushTimer = setTimeout(flushLog, 200);
    }
  };

  eventSource.onerror = (error) =>{
    console.error("Error in streaming output:", error);
    if (gen !== window.__bactflowStreamGen) {
      return;
    }
    flushLog();
    const logText = (document.getElementById("output-bactflow") || {}).innerText || "";
    const hadOutput = /Nextflow|ERROR|Process completed|exited with code|does not exist|N E X T F L O W/i.test(logText);
    if (!hadOutput) {
      BactflowProcessEta.finalizeAll(BactflowTerminal, "stopped");
      BactflowTerminal.append("Stream disconnected.", true);
      BactflowTerminal.setStatus("Stream disconnected", "warn");
    }
    try {
      eventSource.close();
    } catch (e) { /* ignore */ }
    eventSource = null;
    updateButtonStates("stopped");
  };
}

//disconnect function
function disconnectStream(){
  if(eventSource){
    try {
      eventSource.onmessage = null;
      eventSource.onerror = null;
      eventSource.close();
    } catch (e) { /* ignore */ }
    eventSource = null;
  }
  if (!window.__bactflowStreamGen) {
    window.__bactflowStreamGen = 0;
  }
  window.__bactflowStreamGen += 1;
}

function applyUserStopUI() {
  window.__bactflowUserStopped = true;
  disconnectStream();
  BactflowProcessEta.finalizeAll(BactflowTerminal, "stopped");
  updateButtonStates("stopped");
  BactflowTerminal.clear();
  BactflowTerminal.append(STOP_USER_MSG, true);
  BactflowTerminal.setStatus(STOP_USER_MSG, "warn");
}

function clearResultPanels() {
  const panelIds = [
    "quastDiv",
    "baktaDiv",
    "circDiv",
    "linearDiv",
    "circFastaDiv",
    "taxa_class",
    "gtdbTreeDiv",
    "checkmDiv",
    "checkmTreeDiv",
    "plasmidDiv",
    "resultsStatusDiv",
    "rd-results-empty",
    "snpsDiv",
    "svs",
    "abundDiv",
    "abundImgDiv",
    "prevDiv",
    "prevImgDiv",
  ];
  panelIds.forEach((id) => {
    const node = document.getElementById(id);
    if (node) {
      node.style.display = "none";
    }
  });
  const statusEl = document.getElementById("results-status-text");
  if (statusEl) {
    statusEl.textContent = "";
  }
  [
    "output-plasmids",
    "output-viruses",
    "plasmid-summary",
    "virus-summary",
    "plasmid-table-hint",
    "virus-table-hint",
  ].forEach((id) => {
    const node = document.getElementById(id);
    if (node) {
      node.innerHTML = "";
      if (id.endsWith("-hint")) {
        node.textContent = "";
      }
    }
  });
  ["plasmid-section", "virus-section"].forEach((id) => {
    const node = document.getElementById(id);
    if (node) {
      node.style.display = "none";
    }
  });
  [
    "plasmid-chart-by-genome",
    "plasmid-chart-scores",
    "plasmid-chart-lengths",
    "plasmid-chart-topology",
    "virus-chart-by-genome",
    "virus-chart-scores",
    "virus-chart-lengths",
    "virus-chart-topology",
  ].forEach((id) => {
    const node = document.getElementById(id);
    if (node) {
      node.innerHTML = "";
    }
  });
}

function prepareFreshRunUI() {
  window.__bactflowUserStopped = false;
  disconnectStream();
  clearResultPanels();
  BactflowTerminal.clear();
  BactflowTerminal.setStatus("Starting...", "run");
  quastShown = false;
  baktaShown = false;
  taxShown = false;
  checkmShown = false;
  checkmTreeShown = false;
  gtdbTreeShown = false;
  plasmidShown = false;
  resultsAutoShown = false;
  runFinishedSuccess = false;
  setTerminalVisible(true);
  setResultsOpen(false);
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


// Run BactFlow
function run_wf(action){
  if (action === "stop" && window.__bactflowStopInFlight) {
    return;
  }
  if (action === "stop") {
    window.__bactflowStopInFlight = true;
    setTimeout(() => { window.__bactflowStopInFlight = false; }, 1500);
  }

  // including the form data for input parameters to be passed to backend
  const form = document.getElementById('postForm');
  const formData = new FormData(form);

  formData.append("action-assem", action);



  const el = document.getElementById("output-div");
  if(el.style.display ==="none" || el.style.display === ""){
    el.style.display = "block";
  };
  jumpToSection("output-div");
  // switch actions
  switch (action){
    case "run":
      {
        prepareFreshRunUI();
        updateButtonStates("running");
        fetch(`/run_bactflow?action-assem=${action}`, 
          { method: "POST", 
            body : formData }
          )
        .then(async (response) => {
          const message = await response.text();
          if (!response.ok) {
            const looksHtml = /^\s*<!doctype html/i.test(message) || /Werkzeug Debugger/i.test(message);
            BactflowTerminal.clear();
            BactflowTerminal.append(
              looksHtml
                ? "Error starting BactFlow (server error). Restart the post-assembly app and try again."
                : (message || "Error starting BactFlow."),
              true
            );
            BactflowTerminal.setStatus("Failed to start", "error");
            document.getElementById('run-bt').disabled = false;
            document.getElementById('stop-bt').disabled = true;
            document.getElementById('help-bt').disabled = false;
            return;
          }

          BactflowTerminal.append("Bactflow started :)", true);
          BactflowTerminal.setStatus("Running...", "run");
          connectToStream(action);
        })
    
    .catch((error) => {
      BactflowTerminal.append("Failed to start BactFlow: " + error.message, true);
      BactflowTerminal.setStatus("Failed to start", "error");
      document.getElementById('run-bt').disabled = false;
      document.getElementById('stop-bt').disabled = true;
      document.getElementById('help-bt').disabled = false;
    });
        break;
      }
    
    case "stop":
      {
        applyUserStopUI();
        fetch(`/run_bactflow?action-assem=stop`, { method: "POST" })
          .then(async (response) => {
            const message = (await response.text() || "").trim();
            // Keep the dedicated user-stop line; ignore kill/exit noise from the server body.
            if (
              message &&
              message !== STOP_USER_MSG &&
              !/exited with code|Process completed|Stream disconnected|started/i.test(message)
            ) {
              BactflowTerminal.append(message, true);
            }
            BactflowTerminal.setStatus(STOP_USER_MSG, "warn");
            updateButtonStates("stopped");
            refreshCircularPlotButton();
          })
          .catch((error) => {
            BactflowTerminal.append(STOP_USER_MSG + " (" + error.message + ")", true);
            updateButtonStates("stopped");
          });
        break;
      }
        
      case "help":
      {
        BactflowTerminal.clear();
        document.getElementById('help-bt').disabled = true;
    
   

        // hide quast
        const quastDiv = document.getElementById("quastDiv");
        quastDiv.style.display = "none";

      fetch(`/run_bactflow?action-assem=help`, { method: "POST" })
      .then(async (response) => {
        const message = await response.text();
        if(!response.ok) {
          const looksHtml = /^\s*<!doctype html/i.test(message) || /Werkzeug Debugger/i.test(message);
          BactflowTerminal.append(
            looksHtml
              ? "Error showing help (server error). Restart the post-assembly app and try again."
              : (message || "Error showing help for BactFlow."),
            true
          );
          document.getElementById('run-bt').disabled = false;
          document.getElementById('help-bt').disabled = false;
          return;
        };
        BactflowTerminal.clear();
        BactflowTerminal.append("Bactflow's help menu", true);

        // const action = "help";
        updateButtonStates("running");
        connectToStream("help");
        
        
      })
      .catch((error) => {
        BactflowTerminal.append("Failed to give you BactFlow help!" + error.message, true);
        document.getElementById('run-bt').disabled = false;
        document.getElementById('help-bt').disabled = false;
        
      });
        break;
      }
      
    
  }
}

document.addEventListener("DOMContentLoaded", function () {
  const postForm = document.getElementById("postForm");
  if (!postForm) {
    return;
  }

  postForm.addEventListener("submit", (e) => {
    e.preventDefault();
    e.stopPropagation();

    const action =
      (e.submitter && e.submitter.value) ||
      (e.submitter && e.submitter.getAttribute("value")) ||
      null;
    if (!action) {
      return;
    }
    if (action === "run") {
      prepareFreshRunUI();
      document.getElementById("run-bt").disabled = true;
      document.getElementById("stop-bt").disabled = false;
      document.getElementById("help-bt").disabled = true;
    } else if (action !== "stop") {
      window.__bactflowUserStopped = false;
      BactflowTerminal.clear();
      document.getElementById("run-bt").disabled = true;
      document.getElementById("stop-bt").disabled = false;
      document.getElementById("help-bt").disabled = true;
    }

    run_wf(action);
  });

  // Direct clicks — avoid mis-detecting action when submitter is missing.
  const runBt = document.getElementById("run-bt");
  if (runBt) {
    runBt.addEventListener("click", (e) => {
      e.preventDefault();
      e.stopPropagation();
      prepareFreshRunUI();
      document.getElementById("run-bt").disabled = true;
      document.getElementById("stop-bt").disabled = false;
      document.getElementById("help-bt").disabled = true;
      run_wf("run");
    });
  }

  // Direct Stop click — same strategy as assembly, guaranteed even if submitter is missing.
  const stopBt = document.getElementById("stop-bt");
  if (stopBt) {
    stopBt.addEventListener("click", (e) => {
      e.preventDefault();
      e.stopPropagation();
      run_wf("stop");
    });
  }

  const helpBt = document.getElementById("help-bt");
  if (helpBt) {
    helpBt.addEventListener("click", (e) => {
      e.preventDefault();
      e.stopPropagation();
      window.__bactflowUserStopped = false;
      BactflowTerminal.clear();
      document.getElementById("run-bt").disabled = true;
      document.getElementById("stop-bt").disabled = false;
      document.getElementById("help-bt").disabled = true;
      run_wf("help");
    });
  }

  function process_status_hint_from_ui() {
    return document.getElementById("stop-bt")?.disabled ? "stopped" : "running";
  }

  restoreFormData(); //from local storage
  bindStrainFinderFields();

  postForm.addEventListener("input", () => {
    saveFormData();
    updateStrainFinderFromFields();
  });

  const outDirInput = document.getElementById("out_dir");
  if (outDirInput) {
    outDirInput.addEventListener("change", refreshCircularPlotButton);
    outDirInput.addEventListener("blur", refreshCircularPlotButton);
  }

  const baktaSel = document.getElementById("bakta_annot");
  if (baktaSel) {
    baktaSel.addEventListener("change", updateBaktaUI);
  }

  updateBaktaUI();
  updateStrainFinderFromFields();
  bindGeneTypeDropdown();
  bindVizLayout();
  bindLinearPickers();
});

function restoreFormData(){
  let savedData = localStorage.getItem("bactflowFormData");
  if(savedData){
    let formData = JSON.parse(savedData);
    let form = document.getElementById("postForm");
    if (!form) {
      return;
    }

    Array.from(form.elements).forEach(function (element){
      if (element.name && formData[element.name] !== undefined ) {
        if (element.type === "checkbox" || element.type === "radio"){
          element.checked = formData[element.name];
        } else {
          element.value = formData[element.name];
        }
      }
    });
  }
}

function saveFormData(){
  let form = document.getElementById("postForm");
  let formData = {};

  Array.from(form.elements).forEach(function (element){
    if (element.name) {
      if (element.type === "checkbox" || element.type === "radio"){
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

function showChildPanel(childId, visible) {
  const childDiv = document.getElementById(childId);
  if (!childDiv) {
    return;
  }
  childDiv.style.display = visible ? "block" : "none";
}

let baktaReadyPoll = null;

function stopBaktaReadyPoll() {
  if (baktaReadyPoll) {
    clearInterval(baktaReadyPoll);
    baktaReadyPoll = null;
  }
}

function startBaktaReadyPoll() {
  stopBaktaReadyPoll();
  baktaReadyPoll = setInterval(() => {
    const circBtn = document.getElementById("circ-plot-bt");
    if (circBtn && !circBtn.disabled) {
      stopBaktaReadyPoll();
      return;
    }
    refreshCircularPlotButton();
  }, 5000);
}

function setStrainFinderEnabled(ready, message) {
  const abundBtn = document.getElementById("abund-bt");
  const prevBtn = document.getElementById("prev-bt");
  const hint = document.getElementById("strain-finder-hint");
  if (abundBtn) {
    abundBtn.disabled = !ready;
  }
  if (prevBtn) {
    prevBtn.disabled = !ready;
  }
  if (hint) {
    hint.textContent = message || (ready
      ? "Gene and enzyme paths are set. You can run strain finder."
      : "Fill both the enzyme file and the gene annotation directory to activate strain finder.");
  }
}

function strainFinderFieldsFilled() {
  const enzyme = document.querySelector('#postForm [name="enzyme_loc"]');
  const genes = document.querySelector('#postForm [name="gene_files"]');
  const enzymeVal = (enzyme && enzyme.value ? enzyme.value : "").trim();
  const genesVal = (genes && genes.value ? genes.value : "").trim();
  return enzymeVal.length > 0 && genesVal.length > 0;
}

function updateStrainFinderFromFields() {
  const filled = strainFinderFieldsFilled();
  setStrainFinderEnabled(
    filled,
    filled
      ? "Gene and enzyme paths are set. You can run strain finder."
      : "Fill both the enzyme file and the gene annotation directory to activate strain finder."
  );
}

function bindStrainFinderFields() {
  ["enzyme_loc", "gene_files"].forEach((name) => {
    const el = document.querySelector(`#postForm [name="${name}"]`);
    if (!el) {
      return;
    }
    el.addEventListener("input", updateStrainFinderFromFields);
    el.addEventListener("change", updateStrainFinderFromFields);
    el.addEventListener("blur", updateStrainFinderFromFields);
  });
  updateStrainFinderFromFields();
}

function updateBaktaUI() {
  const baktaSel = document.getElementById("bakta_annot");
  const circBtn = document.getElementById("circ-plot-bt");
  const circHint = document.getElementById("circ-plot-hint");
  if (!baktaSel) {
    return;
  }

  const enabled = baktaSel.value === "true";
  showChildPanel("bakta_dbDiv", enabled);
  if (circBtn) {
    circBtn.disabled = true;
  }
  if (circHint) {
    circHint.textContent = "Checking for Bakta annotation files…";
  }
  updateStrainFinderFromFields();
  refreshCircularPlotButton();
  startBaktaReadyPoll();
}

async function refreshCircularPlotButton() {
  const circBtn = document.getElementById("circ-plot-bt");
  const circHint = document.getElementById("circ-plot-hint");

  const form = document.getElementById("postForm");
  const formData = new FormData(form);

  try {
    const res = await fetch("/check-bakta-ready", { method: "POST", body: formData });
    const data = await res.json();
    const plotReady = data.plot_ready === true || data.ready === true;
    if (circBtn) {
      circBtn.disabled = !plotReady;
    }
    if (circHint) {
      const kinds = (data.plot_kinds && data.plot_kinds.length)
        ? data.plot_kinds.join(", ")
        : ".gbk/.gbff/.gff";
      circHint.textContent = plotReady
        ? `Annotation files ready (${data.plot_count || data.gbk_count || 0} file(s): ${kinds}). You can create the plot.`
        : (data.plot_message || data.message || "Run BactFlow with Bakta enabled first.");
    }
    if (plotReady) {
      stopBaktaReadyPoll();
    }
  } catch (error) {
    if (circBtn) {
      circBtn.disabled = true;
    }
    if (circHint) {
      circHint.textContent = "Could not check Bakta annotation files.";
    }
  }
}

async function fetchJsonSafe(url, formData) {
  try {
    const res = await fetch(url, { method: "POST", body: formData });
    const data = await res.json().catch(() => ({}));
    return data && typeof data === "object" ? data : {};
  } catch (_err) {
    return {};
  }
}

function renderCircFastaList(data, opts) {
  const wrap = document.getElementById("circFastaDiv");
  const list = document.getElementById("output-circ-fastas");
  const hint = document.getElementById("circ-fasta-hint");
  if (!wrap || !list) {
    return;
  }
  const files = Array.isArray(data.files) ? data.files : [];
  const circleOn = opts && opts.circleOn;
  const runDone = opts && opts.runDone;
  if (!data.exists || !files.length) {
    if (circleOn && runDone) {
      wrap.style.display = "block";
      if (hint) {
        hint.textContent = "Circulator finished but no FASTA files were found in out_dir/circulated_fasta. Set Resume to False and run again so outputs can be published.";
      }
      list.innerHTML = "";
    } else {
      wrap.style.display = "none";
    }
    return;
  }
  wrap.style.display = "block";
  if (hint) {
    hint.textContent = data.dir
      ? `${files.length} circulated genome(s) in ${data.dir}`
      : `${files.length} circulated genome(s)`;
  }
  list.innerHTML = files.map((name) => `<li><code>${escapeHtml(name)}</code></li>`).join("");
}

let resultsOpen = false;
let resultsAutoShown = false;
let runFinishedSuccess = false;

function runHasFinishedOk() {
  const dash = window.BactflowRunDashboard;
  const eta = window.BactflowProcessEta;
  if (dash && (dash.failed || dash.stopping || dash.frozenRows)) {
    return false;
  }
  if (eta && eta.allDone()) {
    return true;
  }
  return !!(dash && dash.endedAt && dash.startedAt && !dash.failed);
}

function setTerminalVisible(visible) {
  const logs = document.querySelector(".rd-center .rd-logs");
  if (!logs) {
    return;
  }
  logs.classList.toggle("is-hidden", !visible);
  if (!visible) {
    logs.style.height = "";
    return;
  }
  if (window.BactflowRunDashboard && BactflowRunDashboard.alignLogsWithNextflow) {
    requestAnimationFrame(() => BactflowRunDashboard.alignLogsWithNextflow());
  }
}

function showEmptyTerminal() {
  const log = document.getElementById("output-bactflow");
  if (log) {
    log.replaceChildren();
  }
  if (window.BactflowTerminal) {
    BactflowTerminal.setStatus("Run finished", "ok");
  }
}

function isPipelineRunning() {
  const dash = window.BactflowRunDashboard;
  return !!(dash && typeof dash.runState === "function" && dash.runState(dash.procs()) === "running");
}

function anyResultVisible() {
  const pane = document.getElementById("assem-results-pane");
  if (!pane) {
    return false;
  }
  return [...pane.children].some((el) => {
    if (el.id === "rd-results-empty" || el.id === "resultsStatusDiv") {
      return false;
    }
    if (el.style.display === "none") {
      return false;
    }
    return getComputedStyle(el).display !== "none";
  });
}

function setResultsEmptyVisible(visible) {
  const empty = document.getElementById("rd-results-empty");
  if (empty) {
    empty.hidden = !visible;
    empty.style.display = visible ? "block" : "none";
  }
}

function setResultsOpen(open, opts) {
  resultsOpen = !!open;
  const card = document.getElementById("rd-results");
  const btn = document.getElementById("rep-bt");
  if (card) {
    card.classList.toggle("is-shown", resultsOpen);
  }
  if (btn) {
    btn.textContent = resultsOpen ? "Hide results" : "Display results";
  }
  setResultsEmptyVisible(!!(resultsOpen && opts && opts.empty));
  if (!isPipelineRunning()) {
    setTerminalVisible(!resultsOpen);
    if (!resultsOpen && (runFinishedSuccess || runHasFinishedOk())) {
      showEmptyTerminal();
    }
  }
}

function toggleResults() {
  if (resultsOpen) {
    setResultsOpen(false);
    return;
  }
  Promise.resolve(showReport()).then(() => {
    const has = anyResultVisible();
    if (isPipelineRunning()) {
      if (has) {
        setResultsOpen(true);
      }
      return;
    }
    if (runHasFinishedOk()) {
      runFinishedSuccess = true;
    }
    setResultsOpen(true, { empty: !has });
  });
}

function revealResultsOnSuccess(force) {
  const dash = window.BactflowRunDashboard;
  if (resultsAutoShown) {
    return;
  }
  if (!force && !runHasFinishedOk()) {
    return;
  }
  if (dash && (dash.failed || dash.stopping || dash.frozenRows) && !force) {
    return;
  }
  if (isPipelineRunning()) {
    return;
  }
  if (!anyResultVisible()) {
    return;
  }
  resultsAutoShown = true;
  runFinishedSuccess = true;
  setResultsOpen(true);
}

window.toggleResults = toggleResults;
window.revealResultsOnSuccess = revealResultsOnSuccess;
window.showReport = showReport;
window.createGenomePlot = createGenomePlot;

function showReport(){
  const form = document.getElementById("postForm");
  if (!form) {
    return Promise.resolve(false);
  }
  const formData = new FormData(form);
  const outDir = String(formData.get("out_dir") || "").trim();
  if (!outDir) {
    return Promise.resolve(false);
  }

    return Promise.all([
    fetchJsonSafe("/check-quast", formData),
    fetchJsonSafe("/check-bakta-ready", formData),
    fetchJsonSafe("/circular", formData),
    fetchJsonSafe("/taxa-report", formData),
    fetchJsonSafe("/check-circ", formData),
    fetchJsonSafe("/check-checkm", formData),
    fetchJsonSafe("/check-plasmids", formData),
  ])
    .then(([quastData, baktaReady, circPlt, taxData, circFasta, checkmData, plasmidData]) => {
      const quastDiv = document.getElementById("quastDiv");
      const baktaDiv = document.getElementById("baktaDiv");
      const circDiv  = document.getElementById("circDiv");
      const circSpin = document.getElementById("spin-circ");
      const taxDiv = document.getElementById("taxa_class");
      const gtdbTreeDiv = document.getElementById("gtdbTreeDiv");
      const checkmDiv = document.getElementById("checkmDiv");
      const checkmTreeDiv = document.getElementById("checkmTreeDiv");
      const plasmidDiv = document.getElementById("plasmidDiv");

      const circleOn = String(formData.get("circle_genome") || "") === "true";
      const quastOn = String(formData.get("run_quast") || "") === "true";
      const baktaOn = String(formData.get("bakta_annot") || "") === "true";
      const checkmOn = String(formData.get("run_checkm") || "") === "true";
      const runDone = typeof BactflowProcessEta !== "undefined" && BactflowProcessEta.allDone();
      renderCircFastaList(circFasta || {}, { circleOn, runDone });

      const statusEl = document.getElementById("results-status-text");
      const statusWrap = document.getElementById("resultsStatusDiv");
      if (statusEl && statusWrap && (runDone || quastData.exists || (circFasta && circFasta.exists) || baktaReady.plot_ready || taxData.exists || (checkmData && checkmData.exists) || (plasmidData && plasmidData.exists))) {
        const bits = [];
        if (circleOn) {
          bits.push(circFasta && circFasta.exists
            ? `Circulated FASTAs: ${circFasta.count || (circFasta.files || []).length}`
            : "Circulated FASTAs: missing");
        }
        if (quastOn) {
          bits.push(quastData.exists ? "QUAST: ready" : "QUAST: missing");
        }
        if (baktaOn) {
          bits.push((baktaReady.plot_ready || baktaReady.ready) ? "Bakta: ready" : "Bakta: missing");
        }
        if (checkmOn) {
          bits.push(checkmData && checkmData.exists ? "CheckM: ready" : "CheckM: missing");
        }
        if (plasmidData && plasmidData.exists) {
          bits.push("geNomad: ready");
        }
        if (bits.length) {
          statusEl.textContent = bits.join(" · ");
          statusWrap.style.display = "block";
        }
      }

      if (quastData.exists) {
        if (quastDiv) {
          quastDiv.style.display = "block";
        }
        if (!quastShown) {
          quastShown = true;
          quastReport();
        }
      } else if (quastDiv && !quastShown) {
        quastDiv.style.display = "none";
      }

      const baktaPlotReady = baktaReady.plot_ready === true || baktaReady.ready === true;
      if (baktaPlotReady) {
        if (baktaDiv) {
          baktaDiv.style.display = "block";
        }
        if (!baktaShown) {
          baktaShown = true;
          baktaReport();
        }
      } else if (baktaDiv && !baktaShown) {
        baktaDiv.style.display = "none";
      }

      if (circPlt.plot && circDiv) {
        circDiv.style.display = "block";
        if (circSpin) {
          circSpin.style.display = "none";
        }
        showCircularPlot(circPlt.plot);
      }

      if (taxData.exists) {
        if (taxDiv) {
          taxDiv.style.display = "block";
        }
        if (!taxShown) {
          taxShown = true;
          taxReport();
        }
        if (taxData.has_tree && taxData.newick) {
          if (gtdbTreeDiv) {
            gtdbTreeDiv.style.display = "block";
          }
          if (!gtdbTreeShown) {
            gtdbTreeShown = true;
            renderGtdbTree(taxData);
          }
        } else if (gtdbTreeDiv) {
          gtdbTreeDiv.style.display = "none";
        }
      } else {
        if (taxDiv) {
          taxDiv.style.display = "none";
        }
        if (gtdbTreeDiv) {
          gtdbTreeDiv.style.display = "none";
        }
      }

      if (checkmData && checkmData.exists) {
        if (checkmDiv) {
          checkmDiv.style.display = "block";
        }
        if (!checkmShown) {
          checkmShown = true;
          renderCheckmTable(checkmData);
        }
        if (checkmData.has_tree && checkmData.newick) {
          if (checkmTreeDiv) {
            checkmTreeDiv.style.display = "block";
          }
          if (!checkmTreeShown) {
            checkmTreeShown = true;
            renderCheckmTree(checkmData);
          }
        } else if (checkmTreeDiv) {
          checkmTreeDiv.style.display = "none";
        }
      } else {
        if (checkmDiv) {
          checkmDiv.style.display = "none";
        }
        if (checkmTreeDiv) {
          checkmTreeDiv.style.display = "none";
        }
      }

      if (plasmidData && plasmidData.exists) {
        if (plasmidDiv) {
          plasmidDiv.style.display = "block";
        }
        if (!plasmidShown) {
          plasmidShown = true;
          renderPlasmidTable(plasmidData);
        }
      } else if (plasmidDiv) {
        plasmidDiv.style.display = "none";
      }
      revealResultsOnSuccess();
      return anyResultVisible();
    })
    .catch((error) => {
      console.error("Error checking post-assembly reports:", error);
      return false;
    });
}

function renderPlasmidTable(data) {
  if (window.PlasmidViz && typeof window.PlasmidViz.render === "function") {
    window.PlasmidViz.render(data);
    return;
  }
  const wrap = document.getElementById("plasmidDiv");
  const out = document.getElementById("output-plasmids");
  if (!out || !data || !data.exists || !data.plasmid_table) {
    if (wrap) {
      wrap.style.display = "none";
    }
    return;
  }
  if (wrap) {
    wrap.style.display = "block";
  }
  out.innerHTML = data.plasmid_table;
}

async function taxReport() {
  const form = document.getElementById("postForm");
  const formData = new FormData(form);

  try {
    const taxRes = await fetch("/taxa-report", { method: "POST", body: formData });
    const taxDiv = document.getElementById("taxa_class");
    const taxOutputDiv = document.getElementById("output-taxa");
    const abundOutputDiv = document.getElementById("output-taxa-abund");
    const abundTitle = document.getElementById("taxa-abund-title");
    const abundHint = document.getElementById("taxa-abund-hint");
    const taxHint = document.getElementById("taxa-table-hint");

    if (!taxRes.ok) {
      if (taxDiv) {
        taxDiv.style.display = "none";
      }
      const gtdbTreeDiv = document.getElementById("gtdbTreeDiv");
      if (gtdbTreeDiv) {
        gtdbTreeDiv.style.display = "none";
      }
      return;
    }

    const taxData = await taxRes.json();
    if (taxData.exists && taxData.taxa_table && taxOutputDiv) {
      if (taxDiv) {
        taxDiv.style.display = "block";
      }
      if (taxHint) {
        const n = taxData.n_genomes || "";
        taxHint.textContent = n
          ? `GTDB-Tk classification for ${n} genome(s).`
          : "GTDB-Tk classification table.";
      }
      taxOutputDiv.innerHTML = taxData.taxa_table;
      const hasAbund = Boolean(taxData.abund_table && abundOutputDiv);
      if (hasAbund) {
        abundOutputDiv.style.display = "block";
        if (abundTitle) {
          abundTitle.style.display = "block";
        }
        if (abundHint) {
          abundHint.style.display = "block";
          abundHint.textContent = "Counts and percentages of classified genomes at each GTDB rank.";
        }
        abundOutputDiv.innerHTML = taxData.abund_table;
      } else {
        if (abundOutputDiv) {
          abundOutputDiv.style.display = "none";
          abundOutputDiv.innerHTML = "";
        }
        if (abundTitle) {
          abundTitle.style.display = "none";
        }
        if (abundHint) {
          abundHint.style.display = "none";
        }
      }

      setTimeout(() => {
        if (window.$ && $.fn && $.fn.DataTable) {
          if ($.fn.DataTable.isDataTable("#taxa-tab")) {
            $("#taxa-tab").DataTable().destroy();
          }
          $("#taxa-tab").DataTable({
            paging: true,
            pageLength: 10,
            searching: true,
            ordering: true,
            lengthMenu: [[10, 25, 50, -1], [10, 25, 50, "All"]]
          });
          if (hasAbund) {
            if ($.fn.DataTable.isDataTable("#taxa-abund-tab")) {
              $("#taxa-abund-tab").DataTable().destroy();
            }
            $("#taxa-abund-tab").DataTable({
              paging: true,
              pageLength: 10,
              searching: true,
              ordering: true,
              lengthMenu: [[10, 25, 50, -1], [10, 25, 50, "All"]]
            });
          }
        }
      }, 300);
    } else {
      if (taxDiv) {
        taxDiv.style.display = "none";
      }
      const gtdbTreeDiv = document.getElementById("gtdbTreeDiv");
      if (gtdbTreeDiv) {
        gtdbTreeDiv.style.display = "none";
      }
    }
  } catch (error) {
    console.error("Error fetching taxa report:", error);
  }
}

function renderCheckmTable(data) {
  const wrap = document.getElementById("checkmDiv");
  const out = document.getElementById("output-checkm");
  const hint = document.getElementById("checkm-table-hint");
  if (!out) {
    return;
  }
  if (!data || !data.exists || !data.checkm_table) {
    if (wrap) {
      wrap.style.display = "none";
    }
    return;
  }
  if (wrap) {
    wrap.style.display = "block";
  }
  if (hint) {
    const n = data.n_genomes || "";
    hint.textContent = n
      ? `CheckM lineage QA for ${n} genome(s). Completeness, contamination, and marker lineage are from checkm_lineage.txt.`
      : "CheckM lineage QA from checkm_lineage.txt.";
  }
  out.innerHTML = data.checkm_table;
  setTimeout(() => {
    if (!(window.$ && $.fn && $.fn.DataTable)) {
      return;
    }
    if ($.fn.DataTable.isDataTable("#checkm-tab")) {
      $("#checkm-tab").DataTable().destroy();
    }
    $("#checkm-tab").DataTable({
      paging: true,
      pageLength: 10,
      searching: true,
      ordering: true,
      lengthMenu: [[10, 25, 50, -1], [10, 25, 50, "All"]]
    });
  }, 300);
}

function renderCheckmTree(data) {
  const wrap = document.getElementById("checkmTreeDiv");
  const svg = document.getElementById("checkm-tree-svg");
  const hint = document.getElementById("checkm-tree-hint");
  if (!svg || !data || !data.has_tree || !data.newick) {
    if (wrap) {
      wrap.style.display = "none";
    }
    return;
  }
  if (wrap) {
    wrap.style.display = "block";
  }
  if (hint) {
    const src = data.tree_source === "taxon_tree.newick"
      ? "Pruned CheckM taxon tree (taxon_tree.newick). Tips are species names; boxed labels are taxon ranks; edge numbers are branch lengths."
      : "CheckM genome tree. Tips are species names; boxed labels are taxon ranks when present; edge numbers are branch lengths.";
    hint.textContent = src;
  }
  if (typeof CheckmTreeViz === "undefined") {
    svg.textContent = "Tree viewer failed to load.";
    return;
  }
  CheckmTreeViz.render(svg, data.newick, CheckmTreeViz.getLayout(svg) || "rectangular");
  document.querySelectorAll("#checkmTreeDiv .checkm-tree-toggle").forEach((btn) => {
    btn.classList.toggle("active", btn.getAttribute("data-layout") === (CheckmTreeViz.getLayout(svg) || "rectangular"));
  });
}

function setCheckmTreeLayout(layout) {
  const svg = document.getElementById("checkm-tree-svg");
  if (typeof CheckmTreeViz === "undefined") {
    return;
  }
  CheckmTreeViz.setLayout(layout, svg);
  document.querySelectorAll("#checkmTreeDiv .checkm-tree-toggle").forEach((btn) => {
    btn.classList.toggle("active", btn.getAttribute("data-layout") === layout);
  });
}

function downloadCheckmTreePng() {
  if (typeof CheckmTreeViz === "undefined") {
    return;
  }
  CheckmTreeViz.downloadPng("checkm_taxon_tree.png", document.getElementById("checkm-tree-svg"));
}

function renderGtdbTree(data) {
  const wrap = document.getElementById("gtdbTreeDiv");
  const svg = document.getElementById("gtdb-tree-svg");
  const hint = document.getElementById("gtdb-tree-hint");
  if (!svg || !data || !data.has_tree || !data.newick) {
    if (wrap) {
      wrap.style.display = "none";
    }
    return;
  }
  if (wrap) {
    wrap.style.display = "block";
  }
  if (hint) {
    const src = data.tree_source || "classify.tree";
    hint.textContent = `Pruned GTDB-Tk tree (${src}). Tips are your genomes (species names); boxed labels are taxon ranks; edge numbers are branch lengths.`;
  }
  if (typeof CheckmTreeViz === "undefined") {
    svg.textContent = "Tree viewer failed to load.";
    return;
  }
  CheckmTreeViz.render(svg, data.newick, CheckmTreeViz.getLayout(svg) || "rectangular");
  document.querySelectorAll("#gtdbTreeDiv .gtdb-tree-toggle").forEach((btn) => {
    btn.classList.toggle("active", btn.getAttribute("data-layout") === (CheckmTreeViz.getLayout(svg) || "rectangular"));
  });
}

function setGtdbTreeLayout(layout) {
  const svg = document.getElementById("gtdb-tree-svg");
  if (typeof CheckmTreeViz === "undefined") {
    return;
  }
  CheckmTreeViz.setLayout(layout, svg);
  document.querySelectorAll("#gtdbTreeDiv .gtdb-tree-toggle").forEach((btn) => {
    btn.classList.toggle("active", btn.getAttribute("data-layout") === layout);
  });
}

function downloadGtdbTreePng() {
  if (typeof CheckmTreeViz === "undefined") {
    return;
  }
  CheckmTreeViz.downloadPng("gtdbtk_classify_tree.png", document.getElementById("gtdb-tree-svg"));
}

async function quastReport() {
  const form = document.getElementById("postForm");
  const formData = new FormData(form);
  
  try {
   
    let quastResponse = await fetch("/quast-report", { method: "POST", body: formData });
    if (!quastResponse.ok || quastResponse.status === 204) {
      return;
    }
    let quastBlob = await quastResponse.blob();
    if (!quastBlob || quastBlob.size === 0) {
      return;
    }
    let quastUrl = URL.createObjectURL(quastBlob);
    
    let quastOutputDiv = document.getElementById("output-quast");
    quastOutputDiv.innerHTML = `<iframe src="${quastUrl}" style="width: 100%; height: 100%; border: none;"></iframe>`;
   
    let contigOutputDiv = document.getElementById("contig-quast");
    const contigHeading = contigOutputDiv ? contigOutputDiv.previousElementSibling : null;
    try {
      let contigResponse = await fetch("/contig-report", { method: "POST", body: formData });
      if (contigResponse.ok && contigResponse.status !== 204) {
        let contigBlob = await contigResponse.blob();
        if (contigBlob && contigBlob.size > 0 && contigOutputDiv) {
          let contigUrl = URL.createObjectURL(contigBlob);
          contigOutputDiv.style.display = "";
          if (contigHeading && contigHeading.tagName === "H4") {
            contigHeading.style.display = "";
          }
          contigOutputDiv.innerHTML = `<iframe src="${contigUrl}" style="width: 100%; height: 100%; border: none;"></iframe>`;
        }
      } else if (contigOutputDiv) {
        contigOutputDiv.style.display = "none";
        if (contigHeading && contigHeading.tagName === "H4") {
          contigHeading.style.display = "none";
        }
      }
    } catch (_err) {
      if (contigOutputDiv) {
        contigOutputDiv.style.display = "none";
      }
    }
  } catch (error) {
    console.error("Error:", error);
  }
}

async function baktaReport() {
  const form = document.getElementById("postForm");
  const formData = new FormData(form);
  
  try {
    
    let baktaResponse = await fetch("/check-bakta", { method: "POST", body: formData });

   
    if (!baktaResponse.ok) {
      console.error("❌ Bakta report fetch failed:", baktaResponse.status, baktaResponse.statusText);
      throw new Error(`HTTP error! Status: ${baktaResponse.status}`);
    }

    let baktaData = await baktaResponse.json();

    let baktaOutputDiv = document.getElementById("output-bakta");

   
    if (baktaData.exists && baktaData.count_tab) {
  
      baktaOutputDiv.innerHTML = baktaData.count_tab;

      setTimeout(() => {
        if ($.fn.DataTable) {
          console.log("✅ Initializing DataTable...");
          if ($.fn.DataTable.isDataTable("#bakta-tab")) {
            table.DataTable().destroy();
        }
          $("#bakta-tab").DataTable({
            "paging": true,
            "pageLength": 10,
            "searching": true,
            "ordering": true,
            "lengthMenu": [[10, 25, 50, -1], [10, 25, 50, "All"]],
            "responsive": true
          });
        } else {
          console.warn("⚠️ DataTables is not loaded.");
        }
      }, 500);

    } else {
      console.warn("⚠️ No Bakta table found!");
      baktaOutputDiv.innerHTML = `<p>No gene count data available.</p>`;
    }

  } catch (error) {
    console.error("❌ Error fetching Bakta report:", error);
    alert("An error occurred while fetching the Bakta report.");
  }
}

function showSectionError(id, message) {
  const errorBox = document.getElementById(id);
  if (!errorBox) {
    return;
  }
  errorBox.style.display = "block";
  errorBox.textContent = message || "An error occurred.";
}

function hideSectionError(id) {
  const errorBox = document.getElementById(id);
  if (errorBox) {
    errorBox.style.display = "none";
    errorBox.textContent = "";
  }
}

function jumpToSection(id) {
  const el = typeof id === "string" ? document.getElementById(id) : id;
  if (!el) {
    return;
  }
  const pane =
    document.querySelector(".bf-results-col") ||
    el.closest(".bf-results-col") ||
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

async function readJsonError(res, fallback) {
  try {
    const data = await res.json();
    return { data, error: data.error || data.message || fallback };
  } catch (e) {
    return { data: {}, error: fallback };
  }
}

function downloadPlotImage(imgId, filename) {
  const img = document.getElementById(imgId);
  if (!img || !img.src) {
    return;
  }
  const link = document.createElement("a");
  link.href = img.src;
  link.download = filename;
  document.body.appendChild(link);
  link.click();
  link.remove();
}

function downloadCircularPlot() {
  downloadPlotImage("circularImg", "circular_plot.png");
}

function downloadAbundPlot() {
  downloadPlotImage("abundImg", "requested_genes_abundance.jpeg");
}

function downloadPrevPlot() {
  downloadPlotImage("prevImg", "requested_genes_prevalence.jpeg");
}

function geneTypeBoxes() {
  return Array.from(document.querySelectorAll('#geneTypeMenu input[name="gene_type"]'));
}

function selectedGeneTypes() {
  const picked = geneTypeBoxes().filter((el) => el.checked).map((el) => el.value);
  return picked.length ? picked : ["cds"];
}

function syncGeneTypeToggle() {
  const btn = document.getElementById("geneTypeToggle");
  if (btn) {
    btn.textContent = selectedGeneTypes().join(", ");
  }
}

function bindGeneTypeDropdown() {
  const drop = document.getElementById("geneTypeDrop");
  const btn = document.getElementById("geneTypeToggle");
  const menu = document.getElementById("geneTypeMenu");
  if (!drop || !btn || !menu || drop.dataset.bound) {
    return;
  }
  drop.dataset.bound = "1";
  btn.addEventListener("click", (ev) => {
    ev.preventDefault();
    const open = menu.hidden;
    menu.hidden = !open;
    btn.setAttribute("aria-expanded", open ? "true" : "false");
  });
  document.addEventListener("click", (ev) => {
    if (!drop.contains(ev.target)) {
      menu.hidden = true;
      btn.setAttribute("aria-expanded", "false");
    }
  });
  geneTypeBoxes().forEach((box) => box.addEventListener("change", syncGeneTypeToggle));
  syncGeneTypeToggle();
}

function vizLayout() {
  const picked = document.querySelector('input[name="genome_viz_layout"]:checked');
  return picked ? picked.value : "circular";
}

function syncVizLayoutUI() {
  const layout = vizLayout();
  const circ = document.getElementById("circular-opts");
  const lin = document.getElementById("linear-opts");
  const btn = document.getElementById("circ-plot-bt");
  if (circ) {
    circ.style.display = layout === "circular" ? "block" : "none";
  }
  if (lin) {
    lin.style.display = layout === "linear" ? "block" : "none";
  }
  if (btn && !btn.dataset.busy) {
    btn.textContent = layout === "linear" ? "Create linear plot" : "Create circular plot";
  }
}

function bindVizLayout() {
  document.querySelectorAll('input[name="genome_viz_layout"]').forEach((el) => {
    el.addEventListener("change", syncVizLayoutUI);
  });
  syncVizLayoutUI();
}

function fillSelect(sel, values, selected, allLabel) {
  if (!sel) {
    return;
  }
  const current = selected || sel.value;
  sel.innerHTML = "";
  if (allLabel) {
    const opt = document.createElement("option");
    opt.value = "__all__";
    opt.textContent = allLabel;
    sel.appendChild(opt);
  }
  (values || []).forEach((value) => {
    const opt = document.createElement("option");
    opt.value = value;
    opt.textContent = value;
    sel.appendChild(opt);
  });
  if (current && Array.from(sel.options).some((opt) => opt.value === current)) {
    sel.value = current;
  }
}

function bindLinearPickers() {
  const genomeSel = document.getElementById("linear-genome-sel");
  const contigSel = document.getElementById("linear-contig-sel");
  if (genomeSel && !genomeSel.dataset.bound) {
    genomeSel.dataset.bound = "1";
    genomeSel.addEventListener("change", () => createLinearPlot({ keepOpen: true }));
  }
  if (contigSel && !contigSel.dataset.bound) {
    contigSel.dataset.bound = "1";
    contigSel.addEventListener("change", () => createLinearPlot({ keepOpen: true }));
  }
}

function createGenomePlot() {
  if (vizLayout() === "linear") {
    return createLinearPlot();
  }
  return createCircularPlot();
}

const LINEAR_COLORS = {
  cds: ["#3B82F6", "#93C5FD"],
  trna: ["#F59E0B", "#FCD34D"],
  rrna: ["#E11D48", "#FB7185"],
  tmrna: ["#0D9488", "#5EEAD4"],
  ncrna: ["#7C3AED", "#C4B5FD"],
  repeat_region: ["#DB2777", "#F9A8D4"],
  gene: ["#64748B", "#CBD5E1"],
};

function drawLinearCanvas(host, data) {
  if (host._linearMove) {
    window.removeEventListener("mousemove", host._linearMove);
    window.removeEventListener("mouseup", host._linearUp);
  }
  host.innerHTML = "";
  host.style.position = "relative";
  const canvas = document.createElement("canvas");
  canvas.className = "bf-linear-canvas";
  host.appendChild(canvas);
  const tip = document.createElement("div");
  tip.className = "bf-linear-tip";
  tip.hidden = true;
  host.appendChild(tip);
  const features = data.features || [];
  const length = Math.max(1, Number(data.length) || 1);
  const fwdLanes = Math.max(1, Number(data.fwd_lanes) || 1);
  const revLanes = Math.max(1, Number(data.rev_lanes) || 1);
  let view0 = 0;
  let view1 = length;
  let drag = null;

  const colorOf = (feat) => {
    const pair = LINEAR_COLORS[feat.t] || ["#64748B", "#CBD5E1"];
    return feat.strand === "-" ? pair[1] : pair[0];
  };

  const paint = () => {
    const dpr = window.devicePixelRatio || 1;
    const w = Math.max(320, host.clientWidth || 800);
    const h = 440;
    canvas.width = Math.floor(w * dpr);
    canvas.height = Math.floor(h * dpr);
    canvas.style.width = w + "px";
    canvas.style.height = h + "px";
    const ctx = canvas.getContext("2d");
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    const theme = document.documentElement.getAttribute("data-theme");
    ctx.fillStyle = theme === "light" ? "#f8fafc" : "#1e1e1e";
    ctx.fillRect(0, 0, w, h);
    const padL = 18;
    const padR = 18;
    const padT = 32;
    const padB = 40;
    const mid = padT + (h - padT - padB) * (fwdLanes / (fwdLanes + revLanes));
    const span = Math.max(1, view1 - view0);
    const xOf = (bp) => padL + ((bp - view0) / span) * (w - padL - padR);
    ctx.strokeStyle = "rgba(148,163,184,0.75)";
    ctx.lineWidth = 2;
    ctx.beginPath();
    ctx.moveTo(padL, mid);
    ctx.lineTo(w - padR, mid);
    ctx.stroke();
    const laneH = Math.max(6, Math.min(14, (mid - padT - 10) / fwdLanes));
    features.forEach((feat) => {
      if (feat.e < view0 || feat.s > view1) {
        return;
      }
      const x = xOf(Math.max(feat.s, view0));
      const width = Math.max(1.4, xOf(Math.min(feat.e, view1)) - x);
      const lane = feat.lane || 0;
      const y = feat.strand === "-"
        ? mid + 7 + lane * (laneH + 3)
        : mid - 7 - laneH - lane * (laneH + 3);
      ctx.fillStyle = colorOf(feat);
      ctx.fillRect(x, y, width, laneH);
    });
    (data.boundaries || []).forEach((bound) => {
      if (bound.start < view0 || bound.start > view1) {
        return;
      }
      const x = xOf(bound.start);
      ctx.strokeStyle = "rgba(148,163,184,0.35)";
      ctx.setLineDash([3, 4]);
      ctx.beginPath();
      ctx.moveTo(x, padT);
      ctx.lineTo(x, h - padB);
      ctx.stroke();
      ctx.setLineDash([]);
    });
    ctx.fillStyle = theme === "light" ? "#334155" : "#d4d4d4";
    ctx.font = "12px Segoe UI, sans-serif";
    ctx.fillText("Forward", padL, 18);
    ctx.fillText("Reverse", padL, h - 14);
    ctx.textAlign = "right";
    ctx.fillText(data.title || data.length_label || "", w - padR, 18);
    ctx.textAlign = "left";
  };

  const featureAt = (px, py) => {
    const w = host.clientWidth || 800;
    const h = 440;
    const padL = 18;
    const padR = 18;
    const padT = 32;
    const padB = 40;
    const mid = padT + (h - padT - padB) * (fwdLanes / (fwdLanes + revLanes));
    const span = Math.max(1, view1 - view0);
    const laneH = Math.max(6, Math.min(14, (mid - padT - 10) / fwdLanes));
    const bp = view0 + ((px - padL) / Math.max(1, w - padL - padR)) * span;
    for (let i = features.length - 1; i >= 0; i -= 1) {
      const feat = features[i];
      if (bp < feat.s || bp > feat.e) {
        continue;
      }
      const lane = feat.lane || 0;
      const y = feat.strand === "-"
        ? mid + 7 + lane * (laneH + 3)
        : mid - 7 - laneH - lane * (laneH + 3);
      if (py >= y && py <= y + laneH) {
        return feat;
      }
    }
    return null;
  };

  canvas.addEventListener("wheel", (event) => {
    event.preventDefault();
    const rect = canvas.getBoundingClientRect();
    const frac = (event.clientX - rect.left) / rect.width;
    const span = view1 - view0;
    const center = view0 + frac * span;
    const next = Math.min(length, Math.max(Math.min(800, length), span * (event.deltaY > 0 ? 1.25 : 0.8)));
    view0 = Math.max(0, center - next * frac);
    view1 = Math.min(length, view0 + next);
    view0 = Math.max(0, view1 - next);
    paint();
  }, { passive: false });

  canvas.addEventListener("mousedown", (event) => {
    drag = { x: event.clientX, view0, view1 };
  });
  host._linearUp = () => {
    drag = null;
  };
  host._linearMove = (event) => {
    if (!drag) {
      return;
    }
    const rect = canvas.getBoundingClientRect();
    const span = drag.view1 - drag.view0;
    const dx = ((event.clientX - drag.x) / rect.width) * span;
    view0 = Math.max(0, Math.min(length - span, drag.view0 - dx));
    view1 = view0 + span;
    paint();
  };
  window.addEventListener("mouseup", host._linearUp);
  window.addEventListener("mousemove", host._linearMove);
  canvas.addEventListener("mousemove", (event) => {
    if (drag) {
      return;
    }
    const feat = featureAt(event.offsetX, event.offsetY);
    if (!feat) {
      tip.hidden = true;
      return;
    }
    tip.hidden = false;
    tip.textContent = `${feat.n}  ·  ${feat.t}  ·  ${feat.strand === "-" ? "reverse" : "forward"}${feat.p ? "  ·  " + feat.p : ""}`;
    tip.style.left = Math.min(event.offsetX + 14, (host.clientWidth || 400) - 220) + "px";
    tip.style.top = event.offsetY + 16 + "px";
  });
  canvas.addEventListener("mouseleave", () => {
    tip.hidden = true;
  });
  canvas.addEventListener("dblclick", () => {
    view0 = 0;
    view1 = length;
    paint();
  });
  paint();
}

async function createLinearPlot(opts) {
  const form = document.getElementById("postForm");
  const plotBtn = document.getElementById("circ-plot-bt");
  if (!form) {
    return;
  }
  if (plotBtn && plotBtn.disabled && !(opts && opts.keepOpen)) {
    await refreshCircularPlotButton();
    if (plotBtn.disabled) {
      const hint = document.getElementById("circ-plot-hint");
      alert(hint?.textContent || "Bakta annotation results are not available yet.");
      return;
    }
  }
  const formData = new FormData(form);
  const genomeSel = document.getElementById("linear-genome-sel");
  const contigSel = document.getElementById("linear-contig-sel");
  if (genomeSel && genomeSel.value) {
    formData.set("linear_genome", genomeSel.value);
  }
  if (contigSel && contigSel.value) {
    formData.set("linear_contig", contigSel.value);
  }
  const linearDiv = document.getElementById("linearDiv");
  const spin = document.getElementById("spin-linear");
  const typeLabel = document.getElementById("linear-type-label");
  hideSectionError("linear-error");
  if (linearDiv) {
    linearDiv.style.display = "block";
  }
  if (spin) {
    spin.style.display = "flex";
  }
  if (typeLabel) {
    typeLabel.textContent = `(${selectedGeneTypes().join(", ")})`;
  }
  if (!(opts && opts.keepOpen)) {
    setResultsOpen(true);
    jumpToSection("linearDiv");
  }
  if (plotBtn) {
    plotBtn.dataset.busy = "1";
    plotBtn.disabled = true;
    plotBtn.textContent = "Creating plot…";
  }
  try {
    const res = await fetch("/linear", { method: "POST", body: formData });
    const data = await res.json().catch(() => ({}));
    if (spin) {
      spin.style.display = "none";
    }
    if (!res.ok || !data.plot || !(data.features || data.figure)) {
      const reason = data.reason === "no_bakta_dir"
        ? "Bakta annotation output not found. Run BactFlow with Bakta enabled first."
        : (data.error || "Could not create the linear plot.");
      showSectionError("linear-error", reason);
      return;
    }
    hideSectionError("linear-error");
    fillSelect(genomeSel, data.genomes || [], (data.selected || {}).genome);
    fillSelect(contigSel, data.contigs || [], (data.selected || {}).contig, "All contigs");
    const bar = document.getElementById("linear-plot-bar");
    if (bar) {
      bar.style.display = "flex";
    }
    const plotEl = document.getElementById("linearPlot");
    if (plotEl && data.features) {
      drawLinearCanvas(plotEl, data);
    } else if (plotEl && data.figure && window.Plotly) {
      await Plotly.react(plotEl, data.figure.data, data.figure.layout, {
        responsive: true,
        displaylogo: false,
        modeBarButtonsToRemove: ["lasso2d", "select2d"],
      });
    }
    bindLinearPickers();
  } catch (error) {
    if (spin) {
      spin.style.display = "none";
    }
    showSectionError("linear-error", error.message || "An error occurred while creating the linear plot.");
  } finally {
    if (plotBtn) {
      delete plotBtn.dataset.busy;
      plotBtn.disabled = false;
      syncVizLayoutUI();
    }
  }
}

function showCircularPlot(src) {
  const img = document.getElementById("circularImg");
  const bar = document.getElementById("circ-plot-bar");
  const typeLabel = document.getElementById("circ-type-label");
  if (typeLabel) {
    typeLabel.textContent = `(${selectedGeneTypes().join(", ")})`;
  }
  if (img && src && img.getAttribute("src") !== src) {
    img.src = src;
    img.style.display = "block";
  }
  if (bar) {
    bar.style.display = "flex";
  }
}

async function createCircularPlot() {
  const form = document.getElementById("postForm");
  const formData = new FormData(form);
  const circBtn = document.getElementById("circ-plot-bt");

  if (circBtn && circBtn.disabled) {
    await refreshCircularPlotButton();
    if (circBtn.disabled) {
      const hint = document.getElementById("circ-plot-hint");
      alert(hint?.textContent || "Bakta annotation results are not available yet.");
      return;
    }
  }

  formData.set("generate", "true");

  const circDiv = document.getElementById("circDiv");
  const circSpin = document.getElementById("spin-circ");
  const circImg = document.getElementById("circularImg");
  const circBar = document.getElementById("circ-plot-bar");

  hideSectionError("circ-error");
  setResultsOpen(true);
  circDiv.style.display = "block";
  circSpin.style.display = "flex";
  if (circImg) {
    circImg.style.display = "none";
    circImg.removeAttribute("src");
  }
  if (circBar) {
    circBar.style.display = "none";
  }
  jumpToSection("circDiv");
  if (circBtn) {
    circBtn.disabled = true;
    circBtn.textContent = "Creating plot…";
  }

  try {
    const circRes = await fetch("/circular", { method: "POST", body: formData });
    const parsed = await readJsonError(circRes, `HTTP ${circRes.status}`);
    circSpin.style.display = "none";

    if (!circRes.ok || !parsed.data.plot) {
      const reason = parsed.data.reason === "no_bakta_dir"
        ? "Bakta annotation output not found. Run Circulator with Bakta enabled first."
        : (parsed.error || parsed.data.error || "Could not create the circular plot.");
      showSectionError("circ-error", reason);
      return;
    }
    hideSectionError("circ-error");
    showCircularPlot(parsed.data.plot);
  } catch (error) {
    circSpin.style.display = "none";
    console.error("Error creating circular plot:", error);
    showSectionError("circ-error", error.message || "An error occurred while creating the circular plot.");
  } finally {
    if (circBtn) {
      circBtn.disabled = false;
      delete circBtn.dataset.busy;
      circBtn.disabled = false;
      syncVizLayoutUI();
    }
  }
}

async function circReport(){

  const form = document.getElementById("postForm");
  const formData = new FormData(form);
  let circSpin = document.getElementById("spin-circ");
  circSpin.style.display = "block";
  try {

    let circRes = await fetch("/circular", { method: "POST", body: formData });

   
    if (!circRes.ok) {
      console.error("❌ Circular plot fetch failed:", circRes.status, circRes.statusText);
      throw new Error(`HTTP error! Status: ${circRes.status}`);
    }

    let circData = await circRes.json();

    let circImg = document.getElementById("circularImg");
    let circDiv = document.getElementById("output-circ");
   
    if (circData.plot) {
      
      circSpin.style.display = "none";
      showCircularPlot(circData.plot);

    
    } else {
      console.warn("⚠️ No circular plot was found!");
      circDiv.innerHTML = `<p>No circulated plot available.</p>`;
    }

  } catch (error) {
    console.error("❌ Error fetching Bakta report:", error);
    alert("An error occurred while fetching the Bakta report.");
  }
}

let reportCheckInterval = setInterval(() => {
  const runBt = document.getElementById("run-bt");
  if (runBt && runBt.disabled) {
    showReport();
  }
}, 8000);
let quastShown = false;
let baktaShown = false;
let taxShown = false;
let checkmShown = false;
let checkmTreeShown = false;
let gtdbTreeShown = false;
let plasmidShown = false;


bindGeneTypeDropdown();
bindVizLayout();
bindLinearPickers();

// snp finder 
function showSnpError(message) {
  const errorBox = document.getElementById("snp-error");
  if (!errorBox) {
    return;
  }
  errorBox.style.display = "block";
  errorBox.textContent = message || "An error occurred while running SNP finder.";
}

function hideSnpError() {
  const errorBox = document.getElementById("snp-error");
  if (errorBox) {
    errorBox.style.display = "none";
    errorBox.textContent = "";
  }
}

function escapeHtml(value) {
  return String(value ?? "")
    .replace(/&/g, "&amp;")
    .replace(/</g, "&lt;")
    .replace(/>/g, "&gt;")
    .replace(/"/g, "&quot;");
}

function renderChromSummary(chromCounts) {
  const entries = Object.entries(chromCounts || {});
  if (!entries.length) {
    return "";
  }
  const max = Math.max(...entries.map(([, n]) => Number(n) || 0), 1);
  const bars = entries
    .slice(0, 20)
    .map(([chrom, n]) => {
      const width = Math.max(4, Math.round((Number(n) / max) * 100));
      return `<div class="vcf-chrom-row"><span class="vcf-chrom-name">${escapeHtml(chrom)}</span>
        <span class="vcf-chrom-bar"><span style="width:${width}%"></span></span>
        <span class="vcf-chrom-n">${escapeHtml(n)}</span></div>`;
    })
    .join("");
  return `<div class="vcf-chrom-summary"><strong>Variants per contig</strong>${bars}</div>`;
}

function renderVcfGenomeBlock(genome, tableId) {
  const name = escapeHtml(genome.genome || "genome");
  const path = escapeHtml(genome.path || "");
  const n = genome.n_variants ?? 0;
  const shown = genome.n_shown ?? (genome.rows || []).length;
  const truncated = genome.truncated
    ? ` <span class="text-muted">(showing first ${shown} of ${n})</span>`
    : "";

  if (genome.error) {
    return `<div class="vcf-genome-block">
      <h5>${name}</h5>
      <p class="text-danger">${escapeHtml(genome.error)}</p>
      <p class="text-muted small">File: ${path}</p>
    </div>`;
  }

  if (!n) {
    return `<div class="vcf-genome-block">
      <h5>${name}</h5>
      <p>No variants in this VCF (header only).</p>
      <p class="text-muted small">VCF file: <code>${path}</code></p>
    </div>`;
  }

  const columns = genome.columns || ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"];
  const thead = columns.map((c) => `<th>${escapeHtml(c)}</th>`).join("");
  const tbody = (genome.rows || [])
    .map((row) => {
      const cells = columns.map((c) => `<td>${escapeHtml(row[c] ?? "")}</td>`).join("");
      return `<tr>${cells}</tr>`;
    })
    .join("");

  return `<div class="vcf-genome-block">
    <h5>${name} <span class="badge bg-secondary">${n} variant${n === 1 ? "" : "s"}</span>${truncated}</h5>
    <p class="text-muted small mb-2">VCF: <code>${path}</code></p>
    ${renderChromSummary(genome.chrom_counts)}
    <div class="table-responsive">
      <table id="${escapeHtml(tableId)}" class="table table-sm table-striped table-bordered vcf-table" style="width:100%">
        <thead><tr>${thead}</tr></thead>
        <tbody>${tbody}</tbody>
      </table>
    </div>
  </div>`;
}

function renderVariantTables(container, data, idPrefix) {
  if (!container) {
    return;
  }
  const genomes = data.genomes || [];
  if (!genomes.length) {
    const roots = (data.roots || []).filter(Boolean).join(", ");
    const msg = data.message || (
      roots
        ? `SNP / VCF files are in: ${roots}`
        : "No parsed VCF tables were returned. Restart post-assembly and try again."
    );
    container.innerHTML = `<p>${escapeHtml(msg)}</p>`;
    return;
  }

  container.innerHTML = genomes
    .map((g, i) => renderVcfGenomeBlock(g, `${idPrefix}-tab-${i}`))
    .join("<hr>");

  setTimeout(() => {
    if (!window.$ || !$.fn.DataTable) {
      return;
    }
    genomes.forEach((_, i) => {
      const sel = `#${idPrefix}-tab-${i}`;
      if (!document.querySelector(sel)) {
        return;
      }
      if ($.fn.DataTable.isDataTable(sel)) {
        $(sel).DataTable().destroy();
      }
      $(sel).DataTable({
        paging: true,
        pageLength: 25,
        searching: true,
        ordering: true,
        lengthMenu: [[10, 25, 50, 100, -1], [10, 25, 50, 100, "All"]],
        order: [[1, "asc"]],
      });
    });
  }, 200);
}

async function fetchVariantReport(endpoint, form) {
  const fd = new FormData(form);
  const res = await fetch(endpoint, { method: "POST", body: fd });
  let data = {};
  try {
    data = await res.json();
  } catch (e) {
    data = { exists: false, error: "Invalid report response." };
  }
  if (!res.ok) {
    throw new Error(data.error || data.message || `HTTP ${res.status}`);
  }
  return data;
}

function runSNP(){
  let form = document.getElementById("postForm");
  let formData = new FormData(form);
  let snpsDiv = document.getElementById("snpsDiv");
  let snpsOUT = document.getElementById("snp-result");
  let spin = document.getElementById("spin-snp");
  const snpBtn = document.getElementById("snp-bt");

  snpsDiv.style.display = "block";
  hideSnpError();
  if (snpsOUT) {
    snpsOUT.innerHTML = "";
  }
  if (spin) {
    spin.style.display = "flex";
  }
  if (snpBtn) {
    snpBtn.disabled = true;
  }
  jumpToSection("snpsDiv");

  const finishOk = (data) => {
    if (spin) {
      spin.style.display = "none";
    }
    hideSnpError();
    renderVariantTables(snpsOUT, data, "snp");
    jumpToSection("snpsDiv");
  };

  // Prefer existing VCFs (fast). Only re-run the pipeline if none are found.
  fetchVariantReport("/snp-report", form)
    .then((existing) => {
      if (existing.exists && (existing.genomes || []).length) {
        finishOk(existing);
        return null;
      }
      return fetch("/snp-finder", { method: "POST", body: formData }).then(async (res) => {
        let data = {};
        try {
          data = await res.json();
        } catch (e) {
          data = { exists: false, error: "SNP finder returned an invalid response." };
        }
        if (!res.ok) {
          throw new Error(data.error || `HTTP ${res.status}`);
        }
        return data;
      });
    })
    .then(async (data) => {
      if (data == null) {
        return;
      }
      if (!data.exists) {
        if (spin) {
          spin.style.display = "none";
        }
        showSnpError(data.error || "Your SNP files were not created. Check the reference genome path and try again.");
        return;
      }
      if ((data.genomes || []).length) {
        finishOk(data);
        return;
      }
      try {
        const report = await fetchVariantReport("/snp-report", form);
        finishOk(report);
      } catch (err) {
        const outDir = (formData.get("out_dir") || "").toString();
        finishOk({
          roots: outDir ? [`${outDir.replace(/\/$/, "")}/snps`] : [],
          message: data.message || err.message,
          genomes: [],
        });
      }
    })
    .catch((error) => {
      if (spin) {
        spin.style.display = "none";
      }
      console.error("Error fetching SNPs report:", error);
      showSnpError(error.message || "An error occurred while fetching the SNPs report.");
    })
    .finally(() => {
      if (snpBtn) {
        snpBtn.disabled = false;
      }
    });
}

// SVS
// SVS / variant calling
function showSvsError(message) {
  const errorBox = document.getElementById("svs-error");
  if (!errorBox) {
    return;
  }
  errorBox.style.display = "block";
  errorBox.textContent = message || "An error occurred while running variant calling.";
}

function hideSvsError() {
  const errorBox = document.getElementById("svs-error");
  if (errorBox) {
    errorBox.style.display = "none";
    errorBox.textContent = "";
  }
}

function runVCF(){
  let form = document.getElementById("postForm");
  let formData = new FormData(form);
  let svsDiv = document.getElementById("svs");
  let svsOUT = document.getElementById("svs-result");
  let spin = document.getElementById("spin-svs");
  const vcfBtn = document.getElementById("vcf-bt");

  svsDiv.style.display = "block";
  hideSvsError();
  if (svsOUT) {
    svsOUT.innerHTML = "";
  }
  if (spin) {
    spin.style.display = "flex";
  }
  if (vcfBtn) {
    vcfBtn.disabled = true;
  }
  jumpToSection("svs");

  fetch("/svs-finder", {method: "POST", body: formData})
  .then(async (res) => {
    let data = {};
    try {
      data = await res.json();
    } catch (e) {
      data = { exists: false, error: "Variant calling returned an invalid response." };
    }
    if (!res.ok) {
      throw new Error(data.error || `HTTP ${res.status}`);
    }
    return data;
  })
  .then(async (data) => {
    if (spin) {
      spin.style.display = "none";
    }
    if (!data.exists) {
      showSvsError(data.error || "Variant calling files were not created. Check the reference genome path and try again.");
      return;
    }
    hideSvsError();
    if ((data.genomes || []).length) {
      renderVariantTables(svsOUT, data, "svs");
    } else {
      try {
        const report = await fetchVariantReport("/svs-report", form);
        renderVariantTables(svsOUT, report, "svs");
      } catch (err) {
        const outDir = (formData.get("out_dir") || "").toString();
        renderVariantTables(svsOUT, {
          roots: outDir ? [`${outDir.replace(/\/$/, "")}/vcs`] : [],
          message: data.message || err.message,
          genomes: [],
        }, "svs");
      }
    }
    jumpToSection("svs");
  })
  .catch((error) => {
    if (spin) {
      spin.style.display = "none";
    }
    console.error("Error fetching variant calling report:", error);
    showSvsError(error.message || "An error occurred while fetching the variant calling report.");
  })
  .finally(() => {
    if (vcfBtn) {
      vcfBtn.disabled = false;
    }
  });
}

// abund plot 
async function abundRun(){
  const abundBtn = document.getElementById("abund-bt");
  if (abundBtn && abundBtn.disabled) {
    return;
  }

  const form = document.getElementById("postForm");
  const formData = new FormData(form);
  const abundSpin = document.getElementById("spin-abund");
  const abundImgSpin = document.getElementById("spin-abund-img");
  const abundImg = document.getElementById("abundImgDiv");
  const abundImgOut = document.getElementById("abundImg");
  const abundDiv = document.getElementById("output-abund");
  const abDIV = document.getElementById("abundDiv");

  hideSectionError("abund-error");
  hideSectionError("abund-plot-error");
  const abundPlotBar = document.getElementById("abund-plot-bar");
  if (abundPlotBar) {
    abundPlotBar.style.display = "none";
  }
  if (abDIV) {
    abDIV.style.display = "block";
  }
  if (abundImg) {
    abundImg.style.display = "block";
  }
  if (abundSpin) {
    abundSpin.style.display = "flex";
  }
  if (abundImgSpin) {
    abundImgSpin.style.display = "flex";
  }
  if (abundImgOut) {
    abundImgOut.style.display = "none";
  }
  if (abundBtn) {
    abundBtn.disabled = true;
  }
  jumpToSection("abundDiv");

  try {
    const abundRes = await fetch("/abund-run", { method: "POST", body: formData });
    const parsed = await readJsonError(abundRes, `HTTP ${abundRes.status}`);
    const abundData = parsed.data;

    if (abundSpin) {
      abundSpin.style.display = "none";
    }
    if (abundImgSpin) {
      abundImgSpin.style.display = "none";
    }

    if (!abundRes.ok || !abundData.exists) {
      showSectionError("abund-error", parsed.error || abundData.error || "Abundance strain finder failed.");
      return;
    }

    if (abundData.abund_table && abundDiv) {
      abundDiv.innerHTML = abundData.abund_table;
      setTimeout(() => {
        if ($.fn.DataTable) {
          if ($.fn.DataTable.isDataTable("#abund-tab")) {
            $("#abund-tab").DataTable().destroy();
          }
          $("#abund-tab").DataTable({
            "paging": true,
            "pageLength": 10,
            "searching": true,
            "ordering": true,
            "lengthMenu": [[10, 25, 50, -1], [10, 25, 50, "All"]],
            "responsive": true
          });
        }
      }, 500);
    }

    if (abundData.plot_abund && abundImgOut) {
      abundImgOut.src = abundData.plot_abund;
      abundImgOut.style.display = "block";
      if (abundImg) {
        abundImg.style.display = "block";
      }
      if (abundPlotBar) {
        abundPlotBar.style.display = "flex";
      }
      jumpToSection("abundImgDiv");
    } else {
      if (abundPlotBar) {
        abundPlotBar.style.display = "none";
      }
      showSectionError("abund-plot-error", "Abundance table was created but the plot file was not found.");
    }
  } catch (error) {
    if (abundSpin) {
      abundSpin.style.display = "none";
    }
    if (abundImgSpin) {
      abundImgSpin.style.display = "none";
    }
    console.error("Error fetching Abundance report:", error);
    showSectionError("abund-error", error.message || "An error occurred while fetching the Abundance report.");
  } finally {
    updateStrainFinderFromFields();
  }
}

// prev plot 
async function prevRun(){
  const prevBtn = document.getElementById("prev-bt");
  if (prevBtn && prevBtn.disabled) {
    return;
  }

  const form = document.getElementById("postForm");
  const formData = new FormData(form);
  const prevSpin = document.getElementById("spin-prev");
  const prevImg = document.getElementById("prevImgDiv");
  const prevImgOut = document.getElementById("prevImg");
  const prevDiv = document.getElementById("output-prev");
  const prevDIV = document.getElementById("prevDiv");

  hideSectionError("prev-error");
  hideSectionError("prev-plot-error");
  const prevPlotBar = document.getElementById("prev-plot-bar");
  if (prevPlotBar) {
    prevPlotBar.style.display = "none";
  }
  if (prevDIV) {
    prevDIV.style.display = "block";
  }
  if (prevImg) {
    prevImg.style.display = "block";
  }
  if (prevSpin) {
    prevSpin.style.display = "flex";
  }
  if (prevBtn) {
    prevBtn.disabled = true;
  }
  jumpToSection("prevDiv");

  try {
    const prevRes = await fetch("/prev-run", { method: "POST", body: formData });
    const parsed = await readJsonError(prevRes, `HTTP ${prevRes.status}`);
    const prevData = parsed.data;

    if (prevSpin) {
      prevSpin.style.display = "none";
    }

    if (!prevRes.ok || !prevData.exists) {
      showSectionError("prev-error", parsed.error || prevData.error || "Prevalance strain finder failed.");
      return;
    }

    if (prevData.prev_table && prevDiv) {
      prevDiv.innerHTML = prevData.prev_table;
      setTimeout(() => {
        if ($.fn.DataTable) {
          if ($.fn.DataTable.isDataTable("#prev-tab")) {
            $("#prev-tab").DataTable().destroy();
          }
          $("#prev-tab").DataTable({
            "paging": true,
            "pageLength": 10,
            "searching": true,
            "ordering": true,
            "lengthMenu": [[10, 25, 50, -1], [10, 25, 50, "All"]],
            "responsive": true
          });
        }
      }, 500);
    }

    if (prevData.plot_prev && prevImgOut) {
      prevImgOut.src = prevData.plot_prev;
      prevImgOut.style.display = "block";
      if (prevImg) {
        prevImg.style.display = "block";
      }
      if (prevPlotBar) {
        prevPlotBar.style.display = "flex";
      }
      jumpToSection("prevImgDiv");
    } else {
      if (prevPlotBar) {
        prevPlotBar.style.display = "none";
      }
      showSectionError("prev-plot-error", "Prevalance table was created but the plot file was not found.");
    }
  } catch (error) {
    if (prevSpin) {
      prevSpin.style.display = "none";
    }
    console.error("Error fetching prevalance report:", error);
    showSectionError("prev-error", error.message || "An error occurred while fetching the prevalance report.");
  } finally {
    updateStrainFinderFromFields();
  }
}
//reconnect to stream when on assembly


document.addEventListener("DOMContentLoaded", async () =>{
  //  debugger;
  if (window.location.pathname.includes("assembly")){
   const el = document.getElementById("output-div");
   const bactlfowStatus = await fetch("/bactflow_status");
   const data = await bactlfowStatus.json();
          if(data.status === "running"){
              
              if(el.style.display ==="none" || el.style.display === ""){
                el.style.display = "block";
              };
              const runHistory  = await fetch("/bactflow_output");
              const outputData = await runHistory.output.json();
                
                outputDiv.innerText += outputData;
                
               
                updateButtonStates('running');
            
              
          } else if (data.status === "finisehd"){
          const bactFinished = await  fetch("/bactflow_output");
          const outputData = await bactFinished.output;
         
                outputDiv.innerHTML += outputData + "\nBactflow run finished successfully.\n";
                updateButtonStates("stopped");

          } else if(data.status ==="stopped") {
          const bactStopped = await fetch("/bactflow_output");
          const outputData = await bactStopped.output;
          
                                   
              outputDiv.innerHTML += outputData.output + "\n" + "\nBactflow run was stopped.\n";
              updateButtonStates("stopped");
              

          }
    
  }
});

//disconnect
window.addEventListener("beforeunloaded", ()=>{
  disconnectStream();
});


// Handle directory selection
filePicker.addEventListener("click", (event) => {
  const selectedFiles = Array.from(event.target.files);

  if (selectedFiles.length > 0) {
    // Extract the directory path from the first file
    const selectedDirectory = selectedFiles[0].webkitRelativePath.split("/")[0];
    const absolutePath = selectedFiles[0].path || selectedFiles[0].webkitRelativePath.split("/")[0];

    // Update the input field with the absolute directory path
    fastqDirInput.value = absolutePath;
  }
});


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

  //run on loaded page
  toggleDiv();

  motherDiv.addEventListener("change", toggleDiv)
});
};


// slider height 
const sliderW = document.getElementById('width-range');
const spanW = document.getElementById('width-value');
spanW.textContent = sliderW.value;
sliderW.addEventListener("input", ()=>{
  spanW.textContent = sliderW.value;
})


const sliderH = document.getElementById('height-range');
const spanH = document.getElementById('height-value');
spanH.textContent = sliderH.value;
sliderH.addEventListener("input", ()=>{
  spanH.textContent = sliderH.value;
})

//medaka basecaller model toggle
toggler("medaka_polish", "medakaModel");

//coverage parameter toggle
toggler("coverage_filter", "covDiv");

// filter read
toggler(motherId = "nanofilter", childId = "qualDiv");



//gtdbtk
toggler("tax_class", "gtdbtk_dbDiv")

//checkm
toggler("run_checkm", "checkm_dbDiv")

// geNomad plasmids
toggler("run_plasmids", "genomad_dbDiv")
