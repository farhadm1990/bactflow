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
    BactflowProcessEta.finalizeAll(BactflowTerminal, "stopped");
    BactflowTerminal.append("Stream disconnected.", true);
    BactflowTerminal.setStatus("Stream disconnected", "warn");
    eventSource.close();
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
  if (window.BactflowRunDashboard) {
    BactflowRunDashboard.nfFinished = false;
  }
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
            if (window.BactflowRunDashboard && BactflowRunDashboard.markFailed) {
              BactflowRunDashboard.failMessage = (message || "").replace(/\s+/g, " ").slice(0, 180);
              BactflowRunDashboard.markFailed("");
            }
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
  if (!data.exists || !files.length) {
    wrap.style.display = "none";
    list.innerHTML = "";
    if (hint) {
      hint.textContent = "";
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

function isSetupOnlyName(name) {
  const n = String(name || "").toLowerCase().replace(/[^a-z0-9]+/g, "");
  return n === "envsetup" || n === "testify";
}

function lastPipelineProcessDone() {
  const dash = window.BactflowRunDashboard;
  const eta = window.BactflowProcessEta;
  if (dash && (dash.failed || dash.stopping || dash.frozenRows)) {
    return false;
  }
  if (dash && dash.nfFinished && dash.startedAt && !dash.failed) {
    return true;
  }
  if (!eta || !eta.processes || !eta.processes.size) {
    return false;
  }
  const live = [...eta.processes.values()].filter((p) => p.id !== "_executor");
  if (!live.length || live.some((p) => !p.done)) {
    return false;
  }
  if (live.some((p) => !isSetupOnlyName(p.name))) {
    return true;
  }
  const setupOnly = document.getElementById("setup_only");
  return !!(setupOnly && String(setupOnly.value).toLowerCase() === "true");
}

function runHasFinishedOk() {
  const dash = window.BactflowRunDashboard;
  if (dash && (dash.failed || dash.stopping || dash.frozenRows)) {
    return false;
  }
  return lastPipelineProcessDone();
}

function resultStageAllowed(id) {
  const dash = window.BactflowRunDashboard;
  if (!dash || !dash.startedAt) {
    return true;
  }
  if (lastPipelineProcessDone()) {
    return true;
  }
  try {
    const planned = dash.getPlannedStages();
    const rows = dash.fillSkipped(dash.resolveStages(planned, dash.procs()));
    const row = rows.find((r) => r.stage.id === id);
    if (row && row.state === "done" && row.procs.length) {
      return true;
    }
    const allHaveProcs = planned.every((st) => {
      const hit = rows.find((r) => r.stage.id === st.id);
      return hit && hit.procs.length;
    });
    const allSettled = rows.length && rows.every((r) => r.state === "done" || r.state === "fail");
    return !!(allHaveProcs && allSettled);
  } catch (err) {
    return false;
  }
}

function setTerminalVisible(visible) {
  const logs = document.querySelector(".rd-center .rd-logs") || document.querySelector(".rd-logs");
  if (!logs) {
    return;
  }
  logs.classList.toggle("is-hidden", !visible);
  if (!visible) {
    logs.style.height = "";
    logs.style.display = "none";
    return;
  }
  logs.style.display = "";
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
  if (lastPipelineProcessDone()) {
    return false;
  }
  const runBt = document.getElementById("run-bt");
  if (runBt && runBt.disabled) {
    return true;
  }
  const dash = window.BactflowRunDashboard;
  return !!(dash && typeof dash.runState === "function" && dash.runState(dash.procs()) === "running");
}

const RESULT_PANEL_IDS = [
  "circFastaDiv",
  "quastDiv",
  "baktaDiv",
  "circDiv",
  "linearDiv",
  "taxa_class",
  "gtdbTreeDiv",
  "checkmDiv",
  "checkmTreeDiv",
  "plasmidDiv",
  "snpsDiv",
  "abundDiv",
  "abundImgDiv",
  "prevDiv",
  "prevImgDiv",
];

function panelIsShown(el) {
  if (!el) {
    return false;
  }
  if ((el.style && el.style.display) === "none") {
    return false;
  }
  const display = getComputedStyle(el).display;
  return display !== "none" && display !== "";
}

function anyResultVisible() {
  return RESULT_PANEL_IDS.some((id) => panelIsShown(document.getElementById(id)));
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
  if (isPipelineRunning()) {
    setTerminalVisible(true);
    return;
  }
  setTerminalVisible(!resultsOpen);
  if (!resultsOpen && (runFinishedSuccess || runHasFinishedOk())) {
    showEmptyTerminal();
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

let revealInFlight = false;
function revealResultsOnSuccess(force) {
  const dash = window.BactflowRunDashboard;
  if (resultsAutoShown || revealInFlight) {
    return;
  }
  if (dash && (dash.failed || dash.stopping || dash.frozenRows)) {
    return;
  }
  if (isPipelineRunning()) {
    return;
  }
  if (!force && !runHasFinishedOk()) {
    return;
  }
  const openIfReady = () => {
    if (resultsAutoShown || !anyResultVisible()) {
      return false;
    }
    resultsAutoShown = true;
    runFinishedSuccess = true;
    setResultsOpen(true);
    return true;
  };
  if (openIfReady()) {
    return;
  }
  revealInFlight = true;
  Promise.resolve(showReport())
    .then(() => openIfReady())
    .finally(() => {
      revealInFlight = false;
    });
}

window.toggleResults = toggleResults;
window.revealResultsOnSuccess = revealResultsOnSuccess;
window.onPipelineProcessesDone = function () {
  revealResultsOnSuccess(true);
  [1500, 4000, 8000].forEach((ms) => {
    setTimeout(() => {
      if (!resultsAutoShown) {
        revealResultsOnSuccess(true);
      }
    }, ms);
  });
};
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
      const runDone = runHasFinishedOk();
      const allowQuast = resultStageAllowed("quast");
      const allowBakta = resultStageAllowed("bakta");
      const allowCirc = resultStageAllowed("circ");
      const allowGtdb = resultStageAllowed("gtdb");
      const allowCheckm = resultStageAllowed("checkm");
      const allowPlasmid = resultStageAllowed("genomad");
      renderCircFastaList(allowCirc ? (circFasta || {}) : {}, { circleOn, runDone });

      const statusEl = document.getElementById("results-status-text");
      const statusWrap = document.getElementById("resultsStatusDiv");
      if (statusEl && statusWrap && (runDone || (allowQuast && quastData.exists) || (allowCirc && circFasta && circFasta.exists) || (allowBakta && baktaReady.plot_ready) || (allowGtdb && taxData.exists) || (allowCheckm && checkmData && checkmData.exists) || (allowPlasmid && plasmidData && plasmidData.exists))) {
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

      if (allowQuast && quastData.exists) {
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
      if (allowBakta && baktaPlotReady) {
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

      if (allowBakta && circPlt.plot && circDiv) {
        circDiv.style.display = "block";
        if (circSpin) {
          circSpin.style.display = "none";
        }
        showCircularPlot(circPlt.plot);
      }

      if (allowGtdb && taxData.exists) {
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

      if (allowCheckm && checkmData && checkmData.exists) {
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

      if (allowPlasmid && plasmidData && plasmidData.exists) {
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

function linearDownloadName(host, ext) {
  const view = (host && host._linearView) || {};
  const from = Math.max(1, Math.round(view.view0 || 0) + 1);
  const to = Math.max(from, Math.round(view.view1 || view.length || from));
  return `${(host && host._linearName) || "linear_genome"}_${from}-${to}.${ext}`;
}

function downloadLinearPlot() {
  const host = document.getElementById("linearPlot");
  if (!host) {
    return;
  }
  const save = (blob) => {
    if (!blob) {
      return;
    }
    const url = URL.createObjectURL(blob);
    const link = document.createElement("a");
    link.href = url;
    link.download = linearDownloadName(host, "png");
    document.body.appendChild(link);
    link.click();
    link.remove();
    URL.revokeObjectURL(url);
  };
  if (typeof host._linearExportPng === "function") {
    host._linearExportPng(save);
    return;
  }
  const canvas = host._linearCanvas;
  if (canvas) {
    canvas.toBlob(save, "image/png");
  }
}

function downloadLinearHtml() {
  const host = document.getElementById("linearPlot");
  const snap = host && host._linearSnapshot;
  if (!snap) {
    return;
  }
  const theme = document.documentElement.getAttribute("data-theme") || "";
  const runtime = [
    "const LINEAR_COLOR_DEFAULTS = " + JSON.stringify(LINEAR_COLOR_DEFAULTS) + ";",
    "const LINEAR_COLORS = " + JSON.stringify(LINEAR_COLORS) + ";",
    "const LINEAR_TYPE_LABELS = " + JSON.stringify(LINEAR_TYPE_LABELS) + ";",
    "const LINEAR_BASE_COLORS = " + JSON.stringify(LINEAR_BASE_COLORS) + ";",
    "const LINEAR_COMPLEMENT = " + JSON.stringify(LINEAR_COMPLEMENT) + ";",
    foldSearch.toString(),
    linearTypeKey.toString(),
    linearTypeLabel.toString(),
    linearTickStep.toString(),
    linearTickLabel.toString(),
    linearLabelInk.toString(),
    linearColorPair.toString(),
    drawGeneArrow.toString(),
    packLinearLanes.toString(),
    linearTracks.toString(),
    fillLinearLegend.toString(),
    drawLinearCanvas.toString(),
    linearDownloadName.toString(),
    downloadLinearPlot.toString(),
    "const SNAP = " + JSON.stringify(snap) + ";",
    "drawLinearCanvas(document.getElementById('linearPlot'), SNAP.data, SNAP);",
  ].join("\n");
  const html = `<!DOCTYPE html>
<html data-theme="${theme}">
<head>
<meta charset="utf-8">
<title>${(snap.data && snap.data.title) || "Linear genome map"}</title>
<style>
body{margin:16px;font-family:Segoe UI,sans-serif;background:${theme === "light" ? "#f8fafc" : "#111"};color:${theme === "light" ? "#1e293b" : "#ddd"}}
#linearPlot{position:relative;width:100%}
.bf-linear-canvas{width:100%;border-radius:8px;cursor:grab}
.bf-linear-tip{position:absolute;z-index:3;max-width:280px;padding:6px 8px;border-radius:6px;background:rgba(15,23,42,.92);color:#f8fafc;font-size:12px;pointer-events:none}
button{border:0;border-radius:999px;padding:6px 14px;font-weight:600;background:#2dd4bf;color:#111;cursor:pointer}
</style>
</head>
<body>
<p>Scroll to zoom · drag to pan · double-click to reset · <button type="button" onclick="downloadLinearPlot()">Download PNG</button></p>
<div id="linearPlot"></div>
<div id="linear-legend" class="bf-linear-legend"></div>
<script>
${runtime}
</script>
</body>
</html>`;
  const blob = new Blob([html], { type: "text/html" });
  const url = URL.createObjectURL(blob);
  const link = document.createElement("a");
  link.href = url;
  link.download = linearDownloadName(host, "html");
  document.body.appendChild(link);
  link.click();
  link.remove();
  URL.revokeObjectURL(url);
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

function linearPlotIsLive() {
  const div = document.getElementById("linearDiv");
  const host = document.getElementById("linearPlot");
  return !!(div && div.style.display !== "none" && host && host._linearCanvas);
}

function scheduleLinearTypeRefresh() {
  clearTimeout(window.__linearTypeTimer);
  window.__linearTypeTimer = setTimeout(() => {
    if (vizLayout() !== "linear" || !linearPlotIsLive()) {
      return;
    }
    const host = document.getElementById("linearPlot");
    const snap = host && host._linearSnapshot;
    createLinearPlot({ keepOpen: true, quiet: true, restore: snap });
  }, 180);
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
  geneTypeBoxes().forEach((box) => box.addEventListener("change", () => {
    syncGeneTypeToggle();
    scheduleLinearTypeRefresh();
  }));
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
  if (layout === "linear") {
    fillLinearColorEditor();
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
  const picked = Array.isArray(selected) ? selected.filter(Boolean) : (selected ? [selected] : []);
  const previous = sel.multiple
    ? Array.from(sel.selectedOptions || []).map((opt) => opt.value)
    : (sel.value ? [sel.value] : []);
  const want = picked.length ? picked : previous;
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
  const options = Array.from(sel.options);
  if (sel.multiple) {
    let any = false;
    options.forEach((opt) => {
      opt.selected = want.includes(opt.value);
      any = any || opt.selected;
    });
    if (!any && options.length) {
      options[0].selected = true;
    }
  } else if (want[0] && options.some((opt) => opt.value === want[0])) {
    sel.value = want[0];
  }
}

function bindLinearPickers() {
  const genomeSel = document.getElementById("linear-genome-sel");
  const contigSel = document.getElementById("linear-contig-sel");
  const hypo = document.getElementById("linear-hide-hypo");
  if (genomeSel && !genomeSel.dataset.bound) {
    genomeSel.dataset.bound = "1";
    genomeSel.addEventListener("change", () => createLinearPlot({ keepOpen: true }));
  }
  if (contigSel && !contigSel.dataset.bound) {
    contigSel.dataset.bound = "1";
    contigSel.addEventListener("change", () => createLinearPlot({ keepOpen: true }));
  }
  const setupHypo = document.getElementById("linear-hide-hypo-setup");
  const applyHypo = (checked) => {
    if (hypo) {
      hypo.checked = checked;
    }
    if (setupHypo) {
      setupHypo.checked = checked;
    }
    const host = document.getElementById("linearPlot");
    if (host && typeof host._setHideHypo === "function") {
      host._setHideHypo(checked);
    }
  };
  if (hypo && !hypo.dataset.bound) {
    hypo.dataset.bound = "1";
    hypo.addEventListener("change", () => applyHypo(hypo.checked));
  }
  if (setupHypo && !setupHypo.dataset.bound) {
    setupHypo.dataset.bound = "1";
    setupHypo.addEventListener("change", () => applyHypo(setupHypo.checked));
    if (hypo && setupHypo.checked) {
      hypo.checked = true;
    }
  }
  const search = document.getElementById("linear-search");
  if (search && !search.dataset.bound) {
    search.dataset.bound = "1";
    search.addEventListener("input", () => {
      const plot = document.getElementById("linearPlot");
      if (plot && typeof plot._linearSearch === "function") {
        plot._linearSearch(search.value);
      }
    });
    search.addEventListener("keydown", (event) => {
      if (event.key !== "Enter") {
        return;
      }
      event.preventDefault();
      const plot = document.getElementById("linearPlot");
      if (plot && typeof plot._linearSearchNext === "function") {
        plot._linearSearchNext(event.shiftKey ? -1 : 1);
      }
    });
  }
  ["linear-search-prev", "linear-search-next"].forEach((id) => {
    const btn = document.getElementById(id);
    if (!btn || btn.dataset.bound) {
      return;
    }
    btn.dataset.bound = "1";
    btn.addEventListener("click", () => {
      const plot = document.getElementById("linearPlot");
      if (plot && typeof plot._linearSearchNext === "function") {
        plot._linearSearchNext(id.endsWith("prev") ? -1 : 1);
      }
    });
  });
  fillLinearColorEditor();
}

function createGenomePlot() {
  if (vizLayout() === "linear") {
    return createLinearPlot();
  }
  return createCircularPlot();
}

const LINEAR_COLOR_DEFAULTS = {
  cds: ["#E07A5F", "#F2CC8F"],
  trna: ["#E9C46A", "#F6E2A1"],
  rrna: ["#E76F51", "#F4A698"],
  tmrna: ["#2A9D8F", "#8ED3C7"],
  ncrna: ["#6A994E", "#B5D99C"],
  crispr: ["#9B5DE5", "#D4B3FF"],
  "crispr-repeat": ["#C77DFF", "#E0AAFF"],
  repeat_region: ["#4C6EF5", "#A5B4FC"],
  gene: ["#6D6875", "#B8B3BE"],
  orig: ["#F4A261", "#F8D4B0"],
  oric: ["#F4A261", "#F8D4B0"],
};

const LINEAR_COLORS = JSON.parse(JSON.stringify(LINEAR_COLOR_DEFAULTS));

const LINEAR_TYPE_LABELS = {
  cds: "CDS",
  trna: "tRNA",
  rrna: "rRNA",
  tmrna: "tmRNA",
  ncrna: "ncRNA",
  crispr: "CRISPR",
  "crispr-repeat": "CRISPR repeat",
  repeat_region: "Repeat",
  gene: "Gene",
  orig: "oriG",
  oric: "oriC",
};

function linearTypeKey(type) {
  const key = String(type || "").toLowerCase();
  if (key === "ncrna-region") {
    return "ncrna";
  }
  if (key === "crispr-spacer") {
    return "crispr-repeat";
  }
  return key;
}

function loadLinearColors() {
  try {
    const raw = JSON.parse(localStorage.getItem("bactflow_linear_colors") || "{}");
    Object.keys(raw).forEach((key) => {
      if (Array.isArray(raw[key]) && raw[key].length >= 2) {
        LINEAR_COLORS[key] = [raw[key][0], raw[key][1]];
      }
    });
  } catch (err) {
    /* keep defaults */
  }
}

function saveLinearColors() {
  try {
    localStorage.setItem("bactflow_linear_colors", JSON.stringify(LINEAR_COLORS));
  } catch (err) {
    /* ignore */
  }
}

function linearColorPair(type) {
  const key = linearTypeKey(type);
  return LINEAR_COLORS[key] || LINEAR_COLOR_DEFAULTS[key] || ["#6D6875", "#B8B3BE"];
}

function resetLinearColors() {
  Object.keys(LINEAR_COLOR_DEFAULTS).forEach((key) => {
    LINEAR_COLORS[key] = LINEAR_COLOR_DEFAULTS[key].slice();
  });
  saveLinearColors();
}

function fillLinearColorEditor() {
  const host = document.getElementById("linear-color-editor");
  if (!host || host.dataset.ready) {
    return;
  }
  host.dataset.ready = "1";
  const keys = ["cds", "trna", "rrna", "tmrna", "ncrna", "crispr", "crispr-repeat", "repeat_region"];
  host.innerHTML = keys.map((key) => {
    const pair = linearColorPair(key);
    return `<label class="bf-linear-color-row"><span>${linearTypeLabel(key)}</span>
      <input type="color" data-type="${key}" data-strand="+" value="${pair[0]}">
      <input type="color" data-type="${key}" data-strand="-" value="${pair[1]}">
    </label>`;
  }).join("");
  host.addEventListener("input", (event) => {
    const input = event.target;
    if (!input || !input.dataset.type) {
      return;
    }
    const key = input.dataset.type;
    const pair = linearColorPair(key).slice();
    pair[input.dataset.strand === "-" ? 1 : 0] = input.value;
    LINEAR_COLORS[key] = pair;
    saveLinearColors();
    const plot = document.getElementById("linearPlot");
    if (plot && typeof plot._linearPaint === "function") {
      plot._linearPaint();
    }
    fillLinearLegend((plot && plot._linearSnapshot && plot._linearSnapshot.data && plot._linearSnapshot.data.types) || keys);
  });
  const reset = document.getElementById("linear-color-reset");
  if (reset && !reset.dataset.bound) {
    reset.dataset.bound = "1";
    reset.addEventListener("click", () => {
      resetLinearColors();
      host.querySelectorAll("input[type=color]").forEach((input) => {
        const pair = linearColorPair(input.dataset.type);
        input.value = pair[input.dataset.strand === "-" ? 1 : 0];
      });
      const plot = document.getElementById("linearPlot");
      if (plot && typeof plot._linearPaint === "function") {
        plot._linearPaint();
      }
    });
  }
}

loadLinearColors();

const LINEAR_BASE_COLORS = {
  A: "#43AA8B",
  T: "#F94144",
  G: "#F9C74F",
  C: "#577590",
  N: "#94A3B8",
};

const LINEAR_COMPLEMENT = { A: "T", T: "A", G: "C", C: "G", N: "N" };

function linearLabelInk(hex) {
  const raw = String(hex || "").replace("#", "");
  if (raw.length !== 6) {
    return "#111827";
  }
  const r = parseInt(raw.slice(0, 2), 16);
  const g = parseInt(raw.slice(2, 4), 16);
  const b = parseInt(raw.slice(4, 6), 16);
  return (r * 299 + g * 587 + b * 114) / 1000 > 160 ? "#1a1a1a" : "#ffffff";
}

function drawGeneArrow(ctx, x, y, width, height, reverse) {
  const tip = Math.min(height * 0.9, Math.max(4, width * 0.18), Math.max(2, width - 2));
  ctx.beginPath();
  if (!reverse) {
    ctx.moveTo(x, y);
    ctx.lineTo(x + Math.max(0, width - tip), y);
    ctx.lineTo(x + width, y + height / 2);
    ctx.lineTo(x + Math.max(0, width - tip), y + height);
    ctx.lineTo(x, y + height);
  } else {
    ctx.moveTo(x + Math.min(tip, width), y);
    ctx.lineTo(x + width, y);
    ctx.lineTo(x + width, y + height);
    ctx.lineTo(x + Math.min(tip, width), y + height);
    ctx.lineTo(x, y + height / 2);
  }
  ctx.closePath();
  ctx.fill();
  return tip;
}

function packLinearLanes(features) {
  const pack = (list) => {
    const ends = [];
    list.forEach((feat) => {
      let placed = false;
      for (let i = 0; i < ends.length; i += 1) {
        if (feat.s > ends[i] + 8) {
          ends[i] = feat.e;
          feat.lane = i;
          placed = true;
          break;
        }
      }
      if (!placed) {
        feat.lane = ends.length;
        ends.push(feat.e);
      }
    });
    return Math.max(ends.length, 1);
  };
  const fwd = features.filter((feat) => feat.strand !== "-").sort((a, b) => a.s - b.s);
  const rev = features.filter((feat) => feat.strand === "-").sort((a, b) => a.s - b.s);
  return { fwdLanes: pack(fwd), revLanes: pack(rev) };
}

function linearTracks(data) {
  if (data.tracks && data.tracks.length) {
    return data.tracks;
  }
  return [{
    name: data.title || "genome",
    length: data.length,
    length_label: data.length_label,
    fwd_lanes: data.fwd_lanes,
    rev_lanes: data.rev_lanes,
    features: data.features || [],
    boundaries: data.boundaries || [],
    title: data.title,
  }];
}

function selectedLinearGenomes() {
  const sel = document.getElementById("linear-genome-sel");
  if (!sel) {
    return [];
  }
  return Array.from(sel.selectedOptions || []).map((opt) => opt.value).filter(Boolean);
}

function linearTypeLabel(type) {
  const key = linearTypeKey(type);
  return LINEAR_TYPE_LABELS[key] || String(type || "feature");
}

function foldSearch(s) {
  return String(s || "").toLocaleLowerCase();
}

function linearTickStep(span) {
  const raw = Math.max(Number(span) || 1, 1) / 10;
  const pow = Math.pow(10, Math.floor(Math.log10(raw)));
  const norm = raw / pow;
  let step = 10;
  if (norm <= 1) {
    step = 1;
  } else if (norm <= 2) {
    step = 2;
  } else if (norm <= 5) {
    step = 5;
  }
  return Math.max(1, step * pow);
}

function linearTickLabel(value, step) {
  const n = Math.round(value);
  if (step >= 1e6) {
    return (n / 1e6).toFixed(n % 1e6 ? 2 : 0) + " Mb";
  }
  if (step >= 1e3) {
    return (n / 1e3).toFixed(n % 1e3 ? 1 : 0) + " kb";
  }
  return String(n);
}

function fillLinearLegend(types) {
  const el = document.getElementById("linear-legend");
  if (!el) {
    return;
  }
  const seen = [];
  (types || []).forEach((type) => {
    const key = String(type || "").toLowerCase();
    if (key && !seen.includes(key)) {
      seen.push(key);
    }
  });
  if (!seen.length) {
    el.hidden = true;
    el.innerHTML = "";
    return;
  }
  el.hidden = false;
  el.innerHTML = seen.map((type) => {
    const pair = linearColorPair(type);
    return `<span class="bf-linear-legend-item"><i style="background:${pair[0]}"></i><b style="background:${pair[1]}"></b>${linearTypeLabel(type)}</span>`;
  }).join("");
}

function drawLinearCanvas(host, data, restore) {
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

  const tracks = linearTracks(data).map((track) => ({
    ...track,
    features: (track.features || []).map((feat) => ({ ...feat })),
  }));
  const length = Math.max(1, Number(data.length) || Math.max(...tracks.map((t) => Number(t.length) || 1), 1));
  const types = data.types || [...new Set(tracks.flatMap((t) => (t.features || []).map((f) => f.t)))];
  const hypoBox = document.getElementById("linear-hide-hypo");
  let hideHypo = !!(restore && restore.hideHypo) || !!(hypoBox && hypoBox.checked);
  let view0 = restore && restore.view0 != null ? restore.view0 : 0;
  let view1 = restore && restore.view1 != null ? restore.view1 : length;
  let drag = null;
  const seqByTrack = Object.assign({}, (restore && restore.seqByTrack) || data.seqByTrack || {});
  let seqTimer = 0;
  const seqPending = {};
  host._linearCanvas = canvas;
  host._linearName = (data.title || "linear_genome").replace(/[^\w.-]+/g, "_");

  const colorOf = (feat) => {
    const pair = linearColorPair(feat.t);
    return feat.strand === "-" ? pair[1] : pair[0];
  };

  const visibleOf = (track) => (track.features || []).filter((feat) => !(hideHypo && feat.h));

  const featText = (feat, protein) => {
    if (protein && feat.p) {
      return feat.p;
    }
    return feat.n || feat.p || linearTypeLabel(feat.t);
  };

  let searchQuery = "";
  let searchHits = [];
  let searchIdx = 0;
  let aligned = !!(restore && restore.aligned);
  let alignOrigin = Object.assign({}, (restore && restore.alignOrigin) || {});

  const originOf = (track) => {
    if (!aligned) {
      return 0;
    }
    const name = track && track.name ? track.name : track;
    if (!Object.prototype.hasOwnProperty.call(alignOrigin, name)) {
      return 0;
    }
    return Number(alignOrigin[name]) || 0;
  };

  const nucleotideSpan = () => {
    const w = Math.max(320, host.clientWidth || 800);
    const padL = tracks.length > 1 ? 118 : 64;
    return Math.max(36, Math.min(80, Math.floor((w - padL - 18) / 8)));
  };

  const layout = () => {
    const w = Math.max(320, host.clientWidth || 800);
    const padL = tracks.length > 1 ? 118 : 64;
    const padR = 18;
    const padT = 36;
    const padB = 48;
    const span = Math.max(1, view1 - view0);
    const pxPerBp = (w - padL - padR) / span;
    const showBases = pxPerBp >= 4;
    const midBand = showBases ? 64 : 20;
    const laneH = showBases ? 16 : 13;
    const laneGap = 3;
    const header = tracks.length > 1 ? 16 : 0;
    const rows = tracks.map((track) => {
      const packed = packLinearLanes(visibleOf(track));
      const plot = midBand * 2 + packed.fwdLanes * (laneH + laneGap) + packed.revLanes * (laneH + laneGap);
      return { track, packed, plot, header };
    });
    const h = padT + rows.reduce((sum, row) => sum + row.header + row.plot + 18, 0) + padB;
    let y = padT;
    rows.forEach((row) => {
      row.top = y + row.header;
      row.mid = row.top + midBand + row.packed.fwdLanes * (laneH + laneGap);
      y += row.header + row.plot + 18;
    });
    return { w, h, padL, padR, padT, padB, span, pxPerBp, showBases, midBand, laneH, laneGap, rows };
  };

  const xOf = (bp, box, origin) => {
    const org = origin || 0;
    return box.padL + ((bp - org - view0) / box.span) * (box.w - box.padL - box.padR);
  };

  const paintOn = (target, scale, asScreen) => {
    const box = layout();
    const dpr = Math.max(1, Number(scale) || 1);
    target.width = Math.max(1, Math.floor(box.w * dpr));
    target.height = Math.max(1, Math.floor(box.h * dpr));
    if (asScreen) {
      target.style.width = box.w + "px";
      target.style.height = box.h + "px";
    }
    const ctx = target.getContext("2d", { alpha: false });
    if (!ctx) {
      return box;
    }
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.imageSmoothingEnabled = !asScreen;
    const theme = document.documentElement.getAttribute("data-theme");
    const ink = theme === "light" ? "#334155" : "#d4d4d4";
    const muted = theme === "light" ? "#64748B" : "#9ca3af";
    ctx.fillStyle = theme === "light" ? "#f8fafc" : "#1e1e1e";
    ctx.fillRect(0, 0, box.w, box.h);

    ctx.fillStyle = ink;
    ctx.font = "12px Segoe UI, sans-serif";
    ctx.textAlign = "right";
    ctx.fillText(data.title || data.length_label || "", box.w - box.padR, 18);
    ctx.fillStyle = muted;
    ctx.font = "11px ui-monospace, Consolas, monospace";
    ctx.fillText(
      aligned
        ? `aligned  ·  ${Math.max(1, Math.round(view1 - view0))} bp window`
        : `${Math.max(1, Math.round(view0) + 1).toLocaleString()}–${Math.round(view1).toLocaleString()} bp`,
      box.w - box.padR,
      34
    );

    const step = linearTickStep(box.span);
    const minor = step >= 10 ? step / 5 : (step > 1 ? 1 : 0);
    const first = Math.floor(view0 / step) * step;
    const showTickLabels = !box.showBases || step * box.pxPerBp >= 36;
    const minLabelW = 36;

    box.rows.forEach((row) => {
      const trackLen = Math.max(1, Number(row.track.length) || length);
      const origin = originOf(row.track);
      const g0 = origin + view0;
      const g1 = origin + view1;
      ctx.fillStyle = ink;
      ctx.font = "12px Segoe UI, sans-serif";
      ctx.textAlign = "left";
      ctx.textBaseline = "alphabetic";
      if (tracks.length > 1) {
        const label = row.track.name || row.track.title || "genome";
        const hasHit = aligned && Object.prototype.hasOwnProperty.call(alignOrigin, row.track.name);
        ctx.save();
        ctx.beginPath();
        ctx.rect(4, row.mid - 18, box.padL - 10, 32);
        ctx.clip();
        ctx.fillText(label, 8, row.mid - 4);
        ctx.fillStyle = muted;
        ctx.font = "10px ui-monospace, Consolas, monospace";
        ctx.fillText(
          hasHit ? `${(origin + 1).toLocaleString()} bp` : (row.track.length_label || ""),
          8,
          row.mid + 10
        );
        ctx.restore();
      } else {
        ctx.fillText("Forward", 8, row.mid - 16);
        ctx.fillText("Reverse", 8, row.mid + 22);
      }

      visibleOf(row.track).forEach((feat) => {
        const s0 = feat.s - 1;
        if (feat.e < g0 || s0 > g1) {
          return;
        }
        const x = xOf(Math.max(s0, g0), box, origin);
        const width = Math.max(2, xOf(Math.min(feat.e, g1), box, origin) - x);
        const lane = feat.lane || 0;
        const reverse = feat.strand === "-";
        const y = reverse
          ? row.mid + box.midBand + lane * (box.laneH + box.laneGap)
          : row.mid - box.midBand - box.laneH - lane * (box.laneH + box.laneGap);
        const fill = colorOf(feat);
        ctx.fillStyle = fill;
        const tipW = drawGeneArrow(ctx, x, y, width, box.laneH, reverse);
        const isHit = searchQuery && searchHits.some((hit) => hit.feat.s === feat.s && hit.feat.e === feat.e && hit.feat.n === feat.n && hit.track === row.track.name);
        if (isHit) {
          ctx.strokeStyle = "#F9C74F";
          ctx.lineWidth = 2;
          ctx.stroke();
        }
        const text = featText(feat, box.showBases);
        if ((box.showBases || width >= minLabelW) && text) {
          const pad = reverse ? tipW + 3 : 4;
          const maxW = Math.max(4, width - tipW - 7);
          ctx.save();
          ctx.beginPath();
          ctx.rect(x + (reverse ? tipW + 2 : 3), y, maxW, box.laneH);
          ctx.clip();
          ctx.fillStyle = linearLabelInk(fill);
          ctx.font = `600 ${Math.max(8, Math.min(12, box.laneH - 3))}px Segoe UI, sans-serif`;
          ctx.textAlign = "left";
          ctx.textBaseline = "middle";
          ctx.fillText(text, x + pad, y + box.laneH / 2);
          ctx.restore();
        }
      });

      (row.track.boundaries || []).forEach((bound, idx) => {
        if (!idx || bound.start < g0 || bound.start > g1) {
          return;
        }
        const x = xOf(bound.start, box, origin);
        ctx.strokeStyle = "rgba(148,163,184,0.35)";
        ctx.setLineDash([3, 4]);
        ctx.beginPath();
        ctx.moveTo(x, row.top);
        ctx.lineTo(x, row.mid + box.midBand + row.packed.revLanes * (box.laneH + box.laneGap));
        ctx.stroke();
        ctx.setLineDash([]);
      });

      ctx.strokeStyle = muted;
      ctx.fillStyle = ink;
      ctx.font = "11px ui-monospace, Consolas, monospace";
      ctx.textAlign = "center";
      ctx.textBaseline = "top";
      ctx.lineWidth = 1;
      for (let t = first; t <= view1 + step; t += (minor || step)) {
        if (t < view0 || t > view1) {
          continue;
        }
        const genomic = origin + t;
        if (genomic < 0 || genomic > trackLen) {
          continue;
        }
        const x = xOf(genomic, box, origin);
        const major = Math.abs(t / step - Math.round(t / step)) < 1e-6;
        const tick = major ? 7 : 4;
        ctx.beginPath();
        ctx.moveTo(x, row.mid - tick);
        ctx.lineTo(x, row.mid + tick);
        ctx.stroke();
        if (major && showTickLabels) {
          ctx.fillText(linearTickLabel(aligned ? origin + t + 1 : t, step), x, box.showBases ? row.mid + 40 : row.mid + 10);
        }
      }
      ctx.strokeStyle = "rgba(148,163,184,0.9)";
      ctx.lineWidth = 2;
      ctx.beginPath();
      const axis0 = Math.max(0, g0);
      const axis1 = Math.min(trackLen, g1);
      if (axis1 > axis0) {
        ctx.moveTo(xOf(axis0, box, origin), row.mid);
        ctx.lineTo(xOf(axis1, box, origin), row.mid);
        ctx.stroke();
      }

      const seqCache = seqByTrack[row.track.name];
      if (box.showBases && seqCache && seqCache.seq) {
        const from = Math.max(1, Math.floor(g0) + 1);
        const to = Math.min(trackLen, Math.floor(g1));
        ctx.font = `${Math.min(13, Math.max(8, box.pxPerBp * 0.72))}px ui-monospace, Consolas, monospace`;
        ctx.textAlign = "center";
        ctx.textBaseline = "middle";
        for (let bp = from; bp <= to; bp += 1) {
          const i = bp - seqCache.start;
          if (i < 0 || i >= seqCache.seq.length) {
            continue;
          }
          const fwd = seqCache.seq.charAt(i) || "N";
          const rev = LINEAR_COMPLEMENT[fwd] || "N";
          const x = xOf(bp - 0.5, box, origin);
          ctx.fillStyle = LINEAR_BASE_COLORS[fwd] || LINEAR_BASE_COLORS.N;
          ctx.fillText(fwd, x, row.mid - 16);
          ctx.fillStyle = LINEAR_BASE_COLORS[rev] || LINEAR_BASE_COLORS.N;
          ctx.fillText(rev, x, row.mid + 16);
        }
      }
    });

    let lx = box.padL;
    ctx.textAlign = "left";
    ctx.textBaseline = "middle";
    ctx.font = "11px Segoe UI, sans-serif";
    types.forEach((type) => {
      const pair = linearColorPair(type);
      const label = linearTypeLabel(type);
      ctx.fillStyle = pair[0];
      ctx.fillRect(lx, box.h - 28, 10, 10);
      ctx.fillStyle = pair[1];
      ctx.fillRect(lx + 11, box.h - 28, 10, 10);
      ctx.fillStyle = ink;
      ctx.fillText(label, lx + 26, box.h - 23);
      lx += ctx.measureText(label).width + 42;
    });
    ctx.textBaseline = "alphabetic";
    if (asScreen) {
      host._linearView = { view0, view1, length };
      host._linearSnapshot = {
        data,
        view0,
        view1,
        hideHypo,
        seqByTrack,
        aligned,
        alignOrigin,
      };
    }
    return box;
  };

  const paint = () => {
    paintOn(canvas, window.devicePixelRatio || 1, true);
  };

  const exportScaleFor = (box, maxDim) => {
    const w = Math.max(1, box.w);
    const h = Math.max(1, box.h);
    return Math.max(1, Math.min(maxDim / w, maxDim / h, Math.sqrt((maxDim * maxDim) / (w * h))));
  };

  host._linearExportPng = (done) => {
    const box = layout();
    const caps = [16384, 8192, 4096];
    const tryCap = (idx) => {
      const scale = exportScaleFor(box, caps[idx]);
      const bump = document.createElement("canvas");
      bump.width = Math.max(1, Math.floor(box.w * scale));
      bump.height = Math.max(1, Math.floor(box.h * scale));
      if (!bump.getContext("2d") || bump.width < 2 || bump.height < 2) {
        if (idx + 1 < caps.length) {
          tryCap(idx + 1);
          return;
        }
        canvas.toBlob(done, "image/png");
        return;
      }
      paintOn(bump, scale, false);
      bump.toBlob((blob) => {
        if (blob) {
          done(blob);
          return;
        }
        if (idx + 1 < caps.length) {
          tryCap(idx + 1);
          return;
        }
        canvas.toBlob(done, "image/png");
      }, "image/png");
    };
    tryCap(0);
  };

  const featureAt = (px, py) => {
    const box = layout();
    const rel = view0 + ((px - box.padL) / Math.max(1, box.w - box.padL - box.padR)) * box.span;
    for (let r = box.rows.length - 1; r >= 0; r -= 1) {
      const row = box.rows[r];
      const bp = originOf(row.track) + rel;
      const feats = visibleOf(row.track);
      for (let i = feats.length - 1; i >= 0; i -= 1) {
        const feat = feats[i];
        if (bp < feat.s - 1 || bp > feat.e) {
          continue;
        }
        const lane = feat.lane || 0;
        const y = feat.strand === "-"
          ? row.mid + box.midBand + lane * (box.laneH + box.laneGap)
          : row.mid - box.midBand - box.laneH - lane * (box.laneH + box.laneGap);
        if (py >= y && py <= y + box.laneH) {
          return { ...feat, genome: row.track.name };
        }
      }
    }
    return null;
  };

  const requestSeq = () => {
    const box = layout();
    if (box.pxPerBp < 1.2) {
      return;
    }
    const form = document.getElementById("postForm");
    tracks.forEach((track) => {
      const trackLen = Math.max(1, Number(track.length) || length);
      const origin = originOf(track);
      const start = Math.max(1, Math.floor(origin + view0) + 1);
      const end = Math.min(trackLen, Math.ceil(origin + view1));
      if (start > end) {
        return;
      }
      const cached = seqByTrack[track.name];
      if (cached && cached.seq && start >= cached.start && end <= cached.end) {
        return;
      }
      if (!form || seqPending[track.name]) {
        return;
      }
      const pad = Math.min(2400, Math.max(800, (end - start + 1) * 3));
      const reqStart = Math.max(1, start - pad);
      const reqEnd = Math.min(trackLen, end + pad);
      const body = new FormData(form);
      body.set("linear_genome", track.name);
      const contigSel = document.getElementById("linear-contig-sel");
      if (contigSel && contigSel.value) {
        body.set("linear_contig", contigSel.value);
      }
      body.set("start", String(reqStart));
      body.set("end", String(reqEnd));
      seqPending[track.name] = true;
      fetch("/linear-seq", { method: "POST", body })
        .then((res) => (res.ok ? res.json() : Promise.reject(res.status)))
        .then((payload) => {
          if (!payload || payload.seq == null || payload.seq === "") {
            return false;
          }
          seqByTrack[track.name] = {
            start: Number(payload.start) || reqStart,
            end: Number(payload.end) || reqEnd,
            seq: String(payload.seq || "").toUpperCase(),
          };
          paint();
          return true;
        })
        .catch(() => false)
        .then((ok) => {
          seqPending[track.name] = false;
          if (ok) {
            requestSeq();
          }
        });
    });
  };

  const scheduleSeq = (immediate) => {
    clearTimeout(seqTimer);
    if (immediate) {
      requestSeq();
      return;
    }
    seqTimer = setTimeout(requestSeq, 16);
  };

  canvas.addEventListener("wheel", (event) => {
    event.preventDefault();
    const box = layout();
    const rect = canvas.getBoundingClientRect();
    const left = rect.left + rect.width * (box.padL / box.w);
    const plotW = Math.max(1, rect.width * ((box.w - box.padL - box.padR) / box.w));
    const frac = Math.min(1, Math.max(0, (event.clientX - left) / plotW));
    const span = view1 - view0;
    const center = view0 + frac * span;
    const zoomIn = event.deltaY < 0;
    const next = Math.max(1, span * (zoomIn ? (span < 8000 ? 0.62 : 0.72) : 1.25));
    if (aligned) {
      view0 = center - next * frac;
      view1 = view0 + next;
    } else {
      const bounded = Math.min(length, next);
      view0 = Math.max(0, center - bounded * frac);
      view1 = Math.min(length, view0 + bounded);
      view0 = Math.max(0, view1 - bounded);
    }
    paint();
    scheduleSeq(zoomIn);
  }, { passive: false });

  canvas.addEventListener("mousedown", (event) => {
    drag = { x: event.clientX, view0, view1 };
    canvas.style.cursor = "grabbing";
  });
  host._linearUp = () => {
    if (drag) {
      scheduleSeq();
    }
    drag = null;
    canvas.style.cursor = "grab";
  };
  host._linearMove = (event) => {
    if (!drag) {
      return;
    }
    const rect = canvas.getBoundingClientRect();
    const span = drag.view1 - drag.view0;
    const dx = ((event.clientX - drag.x) / rect.width) * span;
    if (aligned) {
      view0 = drag.view0 - dx;
      view1 = view0 + span;
    } else {
      view0 = Math.max(0, Math.min(length - span, drag.view0 - dx));
      view1 = view0 + span;
    }
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
    const genome = feat.genome ? feat.genome + "  ·  " : "";
    tip.textContent = `${genome}${feat.n}  ·  ${linearTypeLabel(feat.t)}  ·  ${feat.strand === "-" ? "reverse" : "forward"}  ·  ${feat.s}-${feat.e}${feat.p ? "  ·  " + feat.p : ""}`;
    tip.style.left = Math.min(event.offsetX + 14, (host.clientWidth || 400) - 220) + "px";
    tip.style.top = event.offsetY + 16 + "px";
  });
  canvas.addEventListener("mouseleave", () => {
    tip.hidden = true;
  });
  canvas.addEventListener("dblclick", () => {
    aligned = false;
    alignOrigin = {};
    view0 = 0;
    view1 = length;
    paint();
  });
  host._linearPaint = paint;
  host._setHideHypo = (value) => {
    hideHypo = !!value;
    paint();
  };
  const updateSearchCount = () => {
    const el = document.getElementById("linear-search-count");
    if (!el) {
      return;
    }
    const genomesHit = new Set(searchHits.map((hit) => hit.track)).size;
    el.textContent = searchHits.length
      ? `${searchIdx + 1}  ·  ${searchHits.length} hits  ·  ${genomesHit} genome${genomesHit === 1 ? "" : "s"}`
      : (searchQuery.trim() ? "0 matches" : "");
  };
  const hitsByTrack = () => {
    const map = {};
    searchHits.forEach((hit) => {
      if (!map[hit.track]) {
        map[hit.track] = [];
      }
      map[hit.track].push(hit);
    });
    Object.keys(map).forEach((name) => {
      map[name].sort((a, b) => a.feat.s - b.feat.s);
    });
    return map;
  };

  const applyAlign = (nth) => {
    const byTrack = hitsByTrack();
    aligned = true;
    alignOrigin = {};
    let maxOcc = 1;
    tracks.forEach((track) => {
      const list = byTrack[track.name] || [];
      maxOcc = Math.max(maxOcc, list.length);
      if (!list.length) {
        return;
      }
      const hit = list[Math.min(Math.max(0, nth), list.length - 1)];
      alignOrigin[track.name] = hit.feat.s - 1;
    });
    const span = nucleotideSpan();
    view0 = -2;
    view1 = span - 2;
    return maxOcc;
  };

  const jumpToHit = () => {
    if (!searchHits.length) {
      return;
    }
    applyAlign(searchIdx);
    paint();
    scheduleSeq();
  };
  host._linearSearch = (query) => {
    searchQuery = String(query || "");
    const q = foldSearch(searchQuery.trim());
    searchHits = [];
    if (q) {
      tracks.forEach((track) => {
        visibleOf(track).forEach((feat) => {
          const hay = foldSearch(`${feat.n} ${feat.p} ${feat.t} ${feat.c || ""} ${track.name}`);
          if (hay.includes(q)) {
            searchHits.push({ feat, track: track.name });
          }
        });
      });
    }
    searchIdx = 0;
    if (!q) {
      aligned = false;
      alignOrigin = {};
      view0 = 0;
      view1 = length;
    }
    updateSearchCount();
    if (searchHits.length) {
      jumpToHit();
    } else {
      paint();
    }
  };
  host._linearSearchNext = (dir) => {
    if (!searchHits.length) {
      host._linearSearch((document.getElementById("linear-search") || {}).value || searchQuery);
      return;
    }
    const byTrack = hitsByTrack();
    const maxOcc = Math.max(1, ...Object.keys(byTrack).map((name) => byTrack[name].length));
    searchIdx = (searchIdx + (dir < 0 ? -1 : 1) + maxOcc) % maxOcc;
    updateSearchCount();
    jumpToHit();
  };
  fillLinearLegend(types);
  paint();
  scheduleSeq();
  const searchEl = document.getElementById("linear-search");
  if (searchEl && searchEl.value) {
    host._linearSearch(searchEl.value);
  }
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
  const pickedGenomes = selectedLinearGenomes();
  if (pickedGenomes.length) {
    formData.set("linear_genome", pickedGenomes.join(","));
  }
  if (contigSel && contigSel.value) {
    formData.set("linear_contig", contigSel.value);
  }
  const quiet = !!(opts && opts.quiet);
  const linearDiv = document.getElementById("linearDiv");
  const spin = document.getElementById("spin-linear");
  const typeLabel = document.getElementById("linear-type-label");
  hideSectionError("linear-error");
  if (linearDiv) {
    linearDiv.style.display = "block";
  }
  if (spin && !quiet) {
    spin.style.display = "flex";
  }
  if (typeLabel) {
    typeLabel.textContent = `(${selectedGeneTypes().join(", ")})`;
  }
  if (!(opts && opts.keepOpen)) {
    setResultsOpen(true);
    jumpToSection("linearDiv");
  }
  if (plotBtn && !quiet) {
    plotBtn.dataset.busy = "1";
    plotBtn.disabled = true;
    plotBtn.textContent = "Creating plot…";
  }
  try {
    const res = await fetch("/linear", { method: "POST", body: formData });
    const data = await res.json().catch(() => ({}));
    if (spin && !quiet) {
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
    fillSelect(genomeSel, data.genomes || [], (data.selected || {}).genomes || (data.selected || {}).genome);
    fillSelect(contigSel, data.contigs || [], (data.selected || {}).contig, "All contigs");
    const bar = document.getElementById("linear-plot-bar");
    if (bar) {
      bar.style.display = "flex";
    }
    const plotEl = document.getElementById("linearPlot");
    bindLinearPickers();
    if (plotEl && (data.tracks || data.features)) {
      const restore = opts && opts.restore;
      drawLinearCanvas(plotEl, data, restore ? {
        view0: restore.view0,
        view1: restore.view1,
        hideHypo: restore.hideHypo,
        seqByTrack: restore.seqByTrack,
        aligned: restore.aligned,
        alignOrigin: restore.alignOrigin,
      } : undefined);
    } else if (plotEl && data.figure && window.Plotly) {
      await Plotly.react(plotEl, data.figure.data, data.figure.layout, {
        responsive: true,
        displaylogo: false,
        modeBarButtonsToRemove: ["lasso2d", "select2d"],
      });
    }
  } catch (error) {
    if (spin && !quiet) {
      spin.style.display = "none";
    }
    showSectionError("linear-error", error.message || "An error occurred while creating the linear plot.");
  } finally {
    if (plotBtn && !quiet) {
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
