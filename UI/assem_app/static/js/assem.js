const outputDiv = document.getElementById("output-bactflow");
const runButton = document.getElementById("run-bt");
const stopButton = document.getElementById("stop-bt");
const helpButton = document.getElementById("help-bt");
let eventSource = null;


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
    // Do not seal stages on a transient SSE drop while the pipeline may still
    // be running — only mark stopped if the server says the job ended.
    const dash = window.BactflowRunDashboard;
    const stillUp = !!(dash && dash.serverRunning);
    if (!stillUp) {
      BactflowProcessEta.finalizeAll(BactflowTerminal, "stopped");
      updateButtonStates("stopped");
    }
    BactflowTerminal.append("Stream disconnected.", true);
    BactflowTerminal.setStatus(stillUp ? "Stream disconnected — pipeline may still be running" : "Stream disconnected", "warn");
    eventSource.close();
    eventSource = null;
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

function clearResultPanels() {
  const panelIds = [
    "quastDiv",
    "baktaDiv",
    "circFastaDiv",
    "taxa_class",
    "gtdbTreeDiv",
    "checkmDiv",
    "checkmTreeDiv",
    "plasmidDiv",
    "resultsStatusDiv",
    "rd-results-empty",
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
  [
    "plasmid-section",
    "virus-section",
  ].forEach((id) => {
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
  disconnectStream();
  clearResultPanels();
  BactflowTerminal.clear();
  if (window.BactflowRunDashboard && typeof BactflowRunDashboard.clearFailState === "function") {
    BactflowRunDashboard.clearFailState();
    BactflowRunDashboard.startedAt = Date.now();
    BactflowRunDashboard.launchAt = BactflowRunDashboard.startedAt;
    BactflowRunDashboard.serverRunning = true;
  }
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

document.getElementById("assemblerDropdown").addEventListener("change", function(){
  updateAssemblerUI();
});

function showEl(el, visible) {
  if (!el) {
    return;
  }
  el.style.display = visible ? "" : "none";
}

function resetFlyeOptions() {
  const coverageFilter = document.getElementById("coverage_filter");
  if (coverageFilter) {
    coverageFilter.value = "false";
    coverageFilter.dispatchEvent(new Event("change"));
  }
  const nanofilter = document.getElementById("nanofilter");
  if (nanofilter) {
    nanofilter.value = "false";
    nanofilter.dispatchEvent(new Event("change"));
  }
  const medakaPolish = document.getElementById("medaka_polish");
  if (medakaPolish) {
    medakaPolish.value = "false";
    medakaPolish.dispatchEvent(new Event("change"));
  }
}

function updateAssemblerUI() {
  const selectedAssembler = document.getElementById("assemblerDropdown").value;
  const assemblers = ["flye", "unicycler", "spades", "pacbio"];
  assemblers.forEach(assembler => {
    const flag = document.getElementById("run_" + assembler);
    if (flag) {
      flag.value = "false";
    }
  });
  if (selectedAssembler !== "none") {
    const selectedFlag = document.getElementById("run_" + selectedAssembler);
    if (selectedFlag) {
      selectedFlag.value = "true";
    }
  }

  const isNone = selectedAssembler === "none";
  const isFlye = selectedAssembler === "flye";
  const isSpades = selectedAssembler === "spades";
  const isUnicycler = selectedAssembler === "unicycler";
  const isPacbio = selectedAssembler === "pacbio";
  const isLongRead = isFlye || isUnicycler || isPacbio;

  const fastqDirGroup = document.getElementById("fastqDirGroup");
  const fastqLabel = document.getElementById("fastqDirLabel");
  const unicyclerDiv = document.getElementById("unicyclerDiv");
  const flyeOntDiv = document.getElementById("flyeOntDiv");
  const pacbioDiv = document.getElementById("pacbioDiv");
  const shortReadInput = document.getElementById("short_read_dir");
  const flyeCoverage = document.getElementById("flyeCoverageOptions");
  const flyePolish = document.getElementById("flyePolishOptions");
  const nanofilterGroup = document.getElementById("nanofilterGroup");
  const circleGenomeGroup = document.getElementById("circleGenomeGroup");
  const concatReadsGroup = document.getElementById("concatReadsGroup");
  const extensionGroup = document.getElementById("extensionGroup");
  const spadesHint = document.getElementById("spadesHint");
  const pacbioHint = document.getElementById("pacbioHint");

  showEl(fastqDirGroup, !isNone);
  showEl(flyeOntDiv, isFlye);
  showEl(unicyclerDiv, isUnicycler);
  showEl(pacbioDiv, isPacbio);
  showEl(spadesHint, isSpades);
  showEl(pacbioHint, isPacbio);
  // Quality/length first in the form; coverage after. Both available for long-read assemblers.
  showEl(nanofilterGroup, isFlye || isUnicycler || isPacbio);
  showEl(flyeCoverage, isFlye || isPacbio);
  showEl(flyePolish, isFlye);
  showEl(circleGenomeGroup, !isNone);
  showEl(concatReadsGroup, isLongRead);
  showEl(extensionGroup, isLongRead);

  if (shortReadInput) {
    shortReadInput.required = isUnicycler;
    if (!isUnicycler) {
      shortReadInput.value = "";
    }
  }

  if (isSpades || isNone) {
    resetFlyeOptions();
  } else if (!isFlye) {
    const medakaPolish = document.getElementById("medaka_polish");
    if (medakaPolish) {
      medakaPolish.value = "false";
      medakaPolish.dispatchEvent(new Event("change"));
    }
  }

  if (isSpades) {
    const concatReads = document.getElementById("concat_reads");
    if (concatReads) {
      concatReads.value = "false";
    }
  }

  if (isPacbio) {
    const concatReads = document.getElementById("concat_reads");
    if (concatReads) {
      concatReads.value = "true";
    }
  }

  if (fastqLabel) {
    if (isUnicycler) {
      fastqLabel.textContent = "Long-read FASTQ directory";
    } else if (isPacbio) {
      fastqLabel.textContent = "PacBio FASTQ directory";
    } else if (isSpades) {
      fastqLabel.textContent = "Illumina FASTQ directory (paired-end)";
    } else {
      fastqLabel.textContent = "FASTQ Directory";
    }
  }
}

function jumpToSection(id) {
  const el = typeof id === "string" ? document.getElementById(id) : id;
  if (!el) {
    return;
  }
  const pane =
    document.getElementById("assem-results-pane") ||
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

function showStartFailure(message) {
  const msg = String(message || "Failed to start")
    .replace(/\s+/g, " ")
    .trim();
  const out = document.getElementById("output-div");
  if (out) {
    out.style.display = "block";
  }
  BactflowTerminal.clear();
  BactflowTerminal.append(msg, true);
  BactflowTerminal.setStatus("Failed to start", "error");
  if (window.BactflowRunDashboard) {
    const dash = window.BactflowRunDashboard;
    dash.failMessage = msg;
    dash.markFailed("");
    if (typeof dash.render === "function") {
      dash.render();
    }
  }
  const runBt = document.getElementById("run-bt");
  const helpBt = document.getElementById("help-bt");
  if (runBt) runBt.disabled = false;
  if (helpBt) helpBt.disabled = false;
  updateButtonStates("stopped");
  jumpToSection("output-div");
}

function clearPreflightFailUi() {
  if (!window.BactflowRunDashboard) {
    return;
  }
  const dash = window.BactflowRunDashboard;
  if (!dash.failed || dash.serverRunning) {
    return;
  }
  if (typeof dash.clearFailState === "function") {
    dash.clearFailState();
  } else {
    dash.failed = false;
    dash.failMessage = "";
    dash.failedProc = "";
    dash.endedAt = null;
    dash.startedAt = null;
    dash.launchAt = null;
    dash.frozenRows = null;
  }
  if (typeof dash.render === "function") {
    dash.render();
  }
  if (window.BactflowTerminal) {
    BactflowTerminal.setStatus("Ready", "idle");
  }
}

async function preflightFastqInputs() {
  const form = document.getElementById("runForm");
  const formData = new FormData(form);
  const assembler = (document.getElementById("assemblerDropdown") || {}).value || "";
  const needsFastq = ["flye", "unicycler", "spades", "pacbio"].includes(assembler);
  if (!needsFastq) {
    return { ok: true };
  }
  const fastqDir = (document.getElementById("fastq_dir") || {}).value || "";
  if (!String(fastqDir).trim()) {
    return { ok: false, message: "Choose a FASTQ directory before starting the run." };
  }
  formData.set("assembler", assembler);
  let response;
  try {
    response = await fetch("/validate-fastq", { method: "POST", body: formData });
  } catch (err) {
    // Network blip — let /run_bactflow enforce the same checks.
    return { ok: true, skipped: true };
  }
  const text = ((await response.text()) || "").trim();
  let data = null;
  if (text) {
    try {
      data = JSON.parse(text);
    } catch (err) {
      data = null;
    }
  }
  // Old server / missing route / HTML error page → do not block with a vague error.
  if (!data && (response.status === 404 || /<!DOCTYPE|<html/i.test(text))) {
    return { ok: true, skipped: true };
  }
  if (!response.ok || !data || !data.ok) {
    const message =
      (data && data.message) ||
      text ||
      "FASTQ directory check failed (no details from server).";
    return { ok: false, message };
  }
  return data;
}

// Run BactFlow
function run_wf(action){

  // including the form data for input parameters to be passed to backend
  const form = document.getElementById('runForm');
  const formData = new FormData(form);

  formData.append("action-assem", action);


  const el = document.getElementById("output-div");
  if(el.style.display ==="none" || el.style.display === ""){
    el.style.display = "block";
  };
  jumpToSection("output-div");
  // switch actions
  switch (action ){
    case "run":
      {
        prepareFreshRunUI();
        updateButtonStates("running");
        fetch(`/run_bactflow?action-assem=${action}`, { method: "POST", body : formData })
        .then(async (response) => {
        if(!response.ok) {
          const message = ((await response.text()) || "").trim();
          showStartFailure(message || "Error starting BactFlow. It might already be running?!");
          return;
        };
        

        BactflowTerminal.append("Bactflow started :)", true);
        BactflowTerminal.setStatus("Running...", "run");
        jumpToSection("output-div");
        connectToStream(action);
      
      
    })
    
    .catch((error) => {
      showStartFailure("Failed to start BactFlow: " + (error && error.message ? error.message : error));
    });
        break;
      }
    
    case "stop":
      {
        jumpToSection("output-div");
        fetch(`/run_bactflow?action-assem=${action}`, {
          method: "POST"
        })
          .then(async (response) => {
            const message = await response.text();
            if (!response.ok) {
              BactflowTerminal.append(message || "Could not stop BactFlow.", true);
              BactflowTerminal.setStatus("Stop failed", "error");
              jumpToSection("output-div");
              return;
            }
            BactflowProcessEta.finalizeAll(BactflowTerminal, "stopped");
            BactflowTerminal.append(message || "Bactflow stopped successfully!", true);
            BactflowTerminal.setStatus("Stopped", "warn");
            disconnectStream();
            updateButtonStates("stopped");
            jumpToSection("output-div");
          })
          .catch((error) => {
            BactflowTerminal.append("Failed to stop BactFlow: " + error.message, true);
            BactflowTerminal.setStatus("Stop failed", "error");
            jumpToSection("output-div");
          });
        break;
      }
        
      case "help":
      {
        revealTerminalForHelp();
        BactflowTerminal.clear();
        document.getElementById('help-bt').disabled = true;
       

        // hide quast
        const quastDiv = document.getElementById("quastDiv");
        quastDiv.style.display = "none";
        jumpToSection("output-div");

      fetch(`/run_bactflow?action-assem=${action}`, { method: "POST" })
      .then((response) => {
        if(!response.ok) {
          BactflowTerminal.append("Error showing help for BactFlow. It might already be running?!", true);
          document.getElementById('run-bt').disabled = false;
          document.getElementById('help-bt').disabled = false;
          jumpToSection("output-div");
          return;
        };
        BactflowTerminal.clear();
        BactflowTerminal.append("Bactflow's help menu", true);
        jumpToSection("output-div");

        // const action = "help";
        updateButtonStates("running");
        connectToStream(action);
        
        
      })
      .catch((error) => {
        BactflowTerminal.append("Failed to give you BactFlow help!" + error.message, true);
        runButton.disabled = false;
        jumpToSection("output-div");
        
      });
        break;
      }
      
    
  }
}

document.getElementById("runForm").addEventListener("submit", async (e) => {
 e.preventDefault();

  const action = e.submitter.value;
  if (action === "run") {
    prepareFreshRunUI();
  } else if (action !== "stop") {
    BactflowTerminal.clear();
  }
  const selectedAssembler = document.getElementById("assemblerDropdown").value;
  const longReads = document.getElementById("fastq_dir").value.trim();
  if (action === "run" && selectedAssembler === "unicycler") {
    const shortReads = document.getElementById("short_read_dir").value.trim();
    if (!longReads || !shortReads) {
      showStartFailure("Unicycler hybrid assembly needs both a long-read path and a short-read path.");
      return;
    }
  }
  if (action === "run" && ["flye", "pacbio", "spades"].includes(selectedAssembler) && !longReads) {
    showStartFailure(
      selectedAssembler === "spades"
        ? "SPAdes needs an Illumina paired-end FASTQ directory (sample_R1 / sample_R2)."
        : "Choose a FASTQ directory before starting the run."
    );
    return;
  }

  if (action === "run") {
    document.getElementById("run-bt").disabled = true;
    document.getElementById("stop-bt").disabled = true;
    document.getElementById("help-bt").disabled = true;
    try {
      const check = await preflightFastqInputs();
      if (!check.ok) {
        showStartFailure(check.message || "FASTQ directory check failed.");
        return;
      }
      if (check.message) {
        BactflowTerminal.append(check.message, true);
      }
    } catch (err) {
      showStartFailure("Could not validate FASTQ path: " + (err && err.message ? err.message : err));
      return;
    }
  }

  if (action !== "stop") {
    document.getElementById('run-bt').disabled = true;
    document.getElementById('stop-bt').disabled = false; 
    document.getElementById('help-bt').disabled = true;
  }

  run_wf(action);
});

// Persist settings across refresh (per-module key so hub apps don't clobber each other).
const FORM_STORAGE_KEY = "bactflowFormData_assem";

function collectFormData(form) {
  const formData = {};
  Array.from(form.elements).forEach((el) => {
    const key = el.name || el.id;
    if (!key || el.type === "file" || el.type === "button" || el.type === "submit") {
      return;
    }
    if (el.type === "radio") {
      if (el.checked) {
        formData[key] = el.value;
      }
      return;
    }
    if (el.type === "checkbox") {
      const group = el.name
        ? form.querySelectorAll(`input[type="checkbox"][name="${CSS.escape(el.name)}"]`)
        : null;
      if (group && group.length > 1) {
        if (!Array.isArray(formData[el.name])) {
          formData[el.name] = [];
        }
        if (el.checked) {
          formData[el.name].push(el.value);
        }
      } else {
        formData[key] = !!el.checked;
      }
      return;
    }
    formData[key] = el.value;
  });
  return formData;
}

function applyFormData(form, formData) {
  Array.from(form.elements).forEach((el) => {
    const key = el.name || el.id;
    if (!key) {
      return;
    }
    if (el.type === "checkbox" && el.name && Array.isArray(formData[el.name])) {
      el.checked = formData[el.name].includes(el.value);
      return;
    }
    if (formData[key] === undefined) {
      return;
    }
    if (el.type === "radio") {
      el.checked = String(el.value) === String(formData[key]);
    } else if (el.type === "checkbox") {
      el.checked = !!formData[key];
    } else {
      el.value = formData[key];
    }
  });
}

function restoreFormData() {
  const form = document.getElementById("runForm");
  if (!form) {
    return;
  }
  let savedData = localStorage.getItem(FORM_STORAGE_KEY);
  if (!savedData) {
    // One-time migrate from the old shared key if it looks like assembly fields.
    const legacy = localStorage.getItem("bactflowFormData");
    if (legacy && legacy.includes("assemblerDropdown")) {
      savedData = legacy;
    }
  }
  if (!savedData) {
    return;
  }
  try {
    applyFormData(form, JSON.parse(savedData));
  } catch (_err) {
    /* ignore corrupt storage */
  }
}

function saveFormData() {
  const form = document.getElementById("runForm");
  if (!form) {
    return;
  }
  try {
    localStorage.setItem(FORM_STORAGE_KEY, JSON.stringify(collectFormData(form)));
  } catch (_err) {
    /* quota / private mode */
  }
}

function clearFormData() {
  localStorage.removeItem(FORM_STORAGE_KEY);
}

document.addEventListener("DOMContentLoaded", function () {
  const form = document.getElementById("runForm");
  if (!form) {
    return;
  }
  restoreFormData();
  updateAssemblerUI();
  // Re-fire selects so show/hide panels match restored values.
  form.querySelectorAll("select").forEach((sel) => {
    sel.dispatchEvent(new Event("change", { bubbles: true }));
  });
  const persist = () => saveFormData();
  form.addEventListener("input", persist);
  form.addEventListener("change", persist);
});

// Silent report lookup: never surface an error if a report is not ready yet.
function hideQuastUi() {
  const quastDiv = document.getElementById("quastDiv");
  if (quastDiv) {
    quastDiv.style.display = "none";
  }
  quastShown = false;
}

function formEl() {
  return document.getElementById("runForm");
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

function escapeHtml(value) {
  return String(value ?? "")
    .replace(/&/g, "&amp;")
    .replace(/</g, "&lt;")
    .replace(/>/g, "&gt;")
    .replace(/"/g, "&quot;");
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
let quastShown = false;
let baktaShown = false;
let taxShown = false;
let checkmShown = false;
let checkmTreeShown = false;
let gtdbTreeShown = false;
let plasmidShown = false;
let reportCheckInterval = null;

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

function revealTerminalForHelp() {
  setTerminalVisible(true);
  const dash = window.BactflowRunDashboard;
  if (dash && typeof dash.scrollTo === "function") {
    requestAnimationFrame(() => dash.scrollTo(".rd-logs"));
    return;
  }
  const logs = document.querySelector(".rd-center .rd-logs") || document.querySelector(".rd-logs");
  if (logs) {
    logs.scrollIntoView({ behavior: "smooth", block: "nearest" });
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
  if (dash && (dash.failed || dash.stopping || dash.frozenRows)) {
    return false;
  }
  if (dash && typeof dash.runState === "function" && dash.runState(dash.procs()) === "running") {
    return true;
  }
  const runBt = document.getElementById("run-bt");
  if (runBt && runBt.disabled && !(dash && (dash.endedAt || dash.nfFinished))) {
    return true;
  }
  return false;
}

const RESULT_PANEL_IDS = [
  "circFastaDiv",
  "quastDiv",
  "baktaDiv",
  "taxa_class",
  "gtdbTreeDiv",
  "checkmDiv",
  "checkmTreeDiv",
  "plasmidDiv",
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
  const running = isPipelineRunning();
  document.documentElement.classList.toggle("rd-live-results", !!(resultsOpen && running));
  if (card) {
    card.classList.toggle("is-shown", resultsOpen);
  }
  if (btn) {
    btn.textContent = resultsOpen ? "Hide results" : "Display results";
  }
  setResultsEmptyVisible(!!(resultsOpen && opts && opts.empty));
  if (running) {
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
      setTerminalVisible(true);
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

function showReport() {
  const form = formEl();
  if (!form) {
    return Promise.resolve(false);
  }
  const outDir = String(new FormData(form).get("out_dir") || "").trim();
  if (!outDir) {
    return Promise.resolve(false);
  }
  // Fresh FormData per request — reusing one body across parallel fetches can drop fields.
  const fd = () => new FormData(form);

  return Promise.all([
    fetchJsonSafe("/check-quast", fd()),
    fetchJsonSafe("/check-bakta-ready", fd()),
    fetchJsonSafe("/taxa-report", fd()),
    fetchJsonSafe("/check-circ", fd()),
    fetchJsonSafe("/check-checkm", fd()),
    fetchJsonSafe("/check-plasmids", fd()),
  ]).then(async ([quastData, baktaReady, taxData, circFasta, checkmData, plasmidData]) => {
    const quastDiv = document.getElementById("quastDiv");
    const baktaDiv = document.getElementById("baktaDiv");
    const taxDiv = document.getElementById("taxa_class");
    const gtdbTreeDiv = document.getElementById("gtdbTreeDiv");
    const checkmDiv = document.getElementById("checkmDiv");
    const checkmTreeDiv = document.getElementById("checkmTreeDiv");
    const plasmidDiv = document.getElementById("plasmidDiv");

    const formData = fd();
    const circleOn = String(formData.get("circle_genome") || "") === "true";
    const quastOn = String(formData.get("run_quast") || "") === "true";
    const baktaOn = String(formData.get("bakta_annot") || "") === "true";
    const checkmOn = String(formData.get("run_checkm") || "") === "true";
    const runDone = typeof BactflowProcessEta !== "undefined" && BactflowProcessEta.allDone();
    if (circFasta && typeof circFasta.exists === "boolean") {
      renderCircFastaList(circFasta, { circleOn, runDone });
    }

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

    if (quastData.exists === true) {
      if (quastDiv) {
        quastDiv.style.display = "block";
      }
      const quastOut = document.getElementById("output-quast");
      const frame = quastOut && quastOut.querySelector("iframe");
      const wantSrc = quastData.report_url || "";
      const hasFrame = !!(frame && frame.getAttribute("src"));
      const stale = !!(frame && wantSrc && frame.getAttribute("src") !== wantSrc);
      if (!hasFrame || stale) {
        await quastReport(quastData);
      }
      quastShown = !!(quastOut && quastOut.querySelector("iframe") && quastOut.querySelector("iframe").getAttribute("src"));
      if (quastDiv) {
        quastDiv.style.display = "block";
      }
    } else if (quastData.exists === false && !quastShown) {
      if (quastDiv) {
        quastDiv.style.display = "none";
      }
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
    } else if (baktaReady.plot_ready === false && baktaReady.ready === false && !baktaShown) {
      if (baktaDiv) {
        baktaDiv.style.display = "none";
      }
    }

    if (taxData.exists === true) {
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
      }
    } else if (taxData.exists === false && !taxShown) {
      if (taxDiv) {
        taxDiv.style.display = "none";
      }
      if (gtdbTreeDiv && !gtdbTreeShown) {
        gtdbTreeDiv.style.display = "none";
      }
    }

    if (checkmData && checkmData.exists === true) {
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
      }
    } else if (checkmData && checkmData.exists === false && !checkmShown) {
      if (checkmDiv) {
        checkmDiv.style.display = "none";
      }
      if (checkmTreeDiv && !checkmTreeShown) {
        checkmTreeDiv.style.display = "none";
      }
    }

    if (plasmidData && plasmidData.exists === true) {
      if (plasmidDiv) {
        plasmidDiv.style.display = "block";
      }
      if (!plasmidShown) {
        plasmidShown = true;
        renderPlasmidTable(plasmidData);
      }
    } else if (plasmidData && plasmidData.exists === false && !plasmidShown) {
      if (plasmidDiv) {
        plasmidDiv.style.display = "none";
      }
    }
    // Only auto-reveal once while a run is finishing — never on idle polls (that blinks the page).
    if (!resultsAutoShown && isPipelineRunning()) {
      revealResultsOnSuccess();
    }
    return anyResultVisible();
  }).catch((error) => {
    console.error("Error checking assembly reports:", error);
    return false;
  });
}

function checkForQuastReport() {
  showReport();
}

function setQuastIframe(hostId, url) {
  const host = document.getElementById(hostId);
  if (!host || !url) {
    return false;
  }
  host.innerHTML = `<iframe src="${url}" style="width: 100%; height: 100%; border: none; background: #fff;"></iframe>`;
  return true;
}

async function quastReport(urls) {
  const form = document.getElementById("runForm");
  const quastDiv = document.getElementById("quastDiv");
  if (!form || !quastDiv) {
    return false;
  }
  const formData = new FormData(form);
  const outDir = String(formData.get("out_dir") || "").trim();
  if (!outDir) {
    return false;
  }

  try {
    let reportUrl = urls && urls.report_url;
    let contigUrl = urls && urls.contig_url;
    if (!reportUrl) {
      const probe = new FormData(form);
      const check = await fetchJsonSafe("/check-quast", probe);
      if (!check.exists) {
        return false;
      }
      reportUrl = check.report_url;
      contigUrl = check.contig_url;
    }
    if (!reportUrl) {
      return false;
    }

    quastDiv.style.display = "block";
    setQuastIframe("output-quast", reportUrl);

    const contigOutputDiv = document.getElementById("contig-quast");
    const contigHeading = contigOutputDiv ? contigOutputDiv.previousElementSibling : null;
    if (contigUrl) {
      // Confirm the contig viewer file exists before embedding.
      try {
        const head = await fetch(contigUrl, { method: "GET" });
        if (head.ok && contigOutputDiv) {
          contigOutputDiv.style.display = "";
          if (contigHeading && contigHeading.tagName === "H4") {
            contigHeading.style.display = "";
          }
          setQuastIframe("contig-quast", contigUrl);
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
    } else if (contigOutputDiv) {
      contigOutputDiv.style.display = "none";
      if (contigHeading && contigHeading.tagName === "H4") {
        contigHeading.style.display = "none";
      }
    }

    return true;
  } catch (err) {
    console.error("quastReport failed:", err);
    return false;
  }
}

async function taxReport() {
  const form = formEl();
  if (!form) {
    return;
  }
  const formData = new FormData(form);
  try {
    const taxRes = await fetch("/taxa-report", { method: "POST", body: formData });
    const taxDiv = document.getElementById("taxa_class");
    const out = document.getElementById("output-taxa");
    const abundOut = document.getElementById("output-taxa-abund");
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
    if (taxData.exists && taxData.taxa_table && out) {
      if (taxDiv) {
        taxDiv.style.display = "block";
      }
      if (taxHint) {
        const n = taxData.n_genomes || "";
        taxHint.textContent = n
          ? `GTDB-Tk classification for ${n} genome(s).`
          : "GTDB-Tk classification table.";
      }
      out.innerHTML = taxData.taxa_table;
      const hasAbund = Boolean(taxData.abund_table && abundOut);
      if (hasAbund) {
        abundOut.style.display = "block";
        if (abundTitle) {
          abundTitle.style.display = "block";
        }
        if (abundHint) {
          abundHint.style.display = "block";
          abundHint.textContent = "Counts and percentages of classified genomes at each GTDB rank.";
        }
        abundOut.innerHTML = taxData.abund_table;
      } else {
        if (abundOut) {
          abundOut.style.display = "none";
          abundOut.innerHTML = "";
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

async function baktaReport() {
  const form = formEl();
  if (!form) {
    return;
  }
  const formData = new FormData(form);
  try {
    const baktaResponse = await fetch("/check-bakta", { method: "POST", body: formData });
    const baktaOutputDiv = document.getElementById("output-bakta");
    if (!baktaResponse.ok || !baktaOutputDiv) {
      return;
    }
    const baktaData = await baktaResponse.json();
    if (baktaData.exists && baktaData.count_tab) {
      baktaOutputDiv.innerHTML = baktaData.count_tab;
      setTimeout(() => {
        if (window.$ && $.fn && $.fn.DataTable) {
          if ($.fn.DataTable.isDataTable("#bakta-tab")) {
            $("#bakta-tab").DataTable().destroy();
          }
          $("#bakta-tab").DataTable({
            paging: true,
            pageLength: 10,
            searching: true,
            ordering: true,
            lengthMenu: [[10, 25, 50, -1], [10, 25, 50, "All"]]
          });
        }
      }, 300);
    } else {
      baktaOutputDiv.innerHTML = "<p>No gene count data available yet.</p>";
    }
  } catch (error) {
    console.error("Error fetching Bakta report:", error);
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
    hint.textContent = data.tree_source === "taxon_tree.newick"
      ? "Pruned CheckM taxon tree (taxon_tree.newick). Tips are species names; boxed labels are taxon ranks; edge numbers are branch lengths."
      : "CheckM genome tree. Tips are species names; boxed labels are taxon ranks when present; edge numbers are branch lengths.";
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

function selectedGeneTypes() {
  const sel = document.getElementById("geneType");
  if (!sel) {
    return ["cds"];
  }
  const picked = Array.from(sel.selectedOptions).map((opt) => opt.value).filter(Boolean);
  return picked.length ? picked : ["cds"];
}

function updateBaktaUI() {
  const baktaSel = document.getElementById("bakta_annot");
  if (!baktaSel) {
    return;
  }
}

// Poll reports only while BactFlow is running — idle polling rewrote the DOM and blinked the page.
reportCheckInterval = setInterval(() => {
  const runBt = document.getElementById("run-bt");
  if (runBt && runBt.disabled && isPipelineRunning()) {
    showReport();
  }
}, 8000);




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


// Handle directory selection for FASTQ input
if (fastqDirInput) {
  fastqDirInput.addEventListener("input", () => {
    clearPreflightFailUi();
  });
  fastqDirInput.addEventListener("change", (event) => {
    clearPreflightFailUi();
    const selectedFiles = Array.from(event.target.files || []);

    if (selectedFiles.length > 0) {
      const absolutePath = selectedFiles[0].path || selectedFiles[0].webkitRelativePath.split("/")[0];
      fastqDirInput.value = absolutePath;
    }
  });
}
const shortReadDirInput = document.getElementById("short_read_dir");
if (shortReadDirInput) {
  shortReadDirInput.addEventListener("input", () => {
    clearPreflightFailUi();
  });
  shortReadDirInput.addEventListener("change", () => {
    clearPreflightFailUi();
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

  //run on loaded page
  toggleDiv();

  motherDiv.addEventListener("change", toggleDiv)
});
};

//medaka basecaller model toggle
toggler("medaka_polish", "medakaModel");

//coverage parameter toggle
toggler("coverage_filter", "covDiv");

// filter read
toggler(motherId = "nanofilter", childId = "qualDiv");

// medaka tensor
toggler(motherId = "medaka_polish", childID = "tesnDiv")

//bakta
toggler("bakta_annot", "bakta_dbDiv")
document.addEventListener("DOMContentLoaded", () => {
  const baktaSel = document.getElementById("bakta_annot");
  if (baktaSel) {
    baktaSel.addEventListener("change", updateBaktaUI);
    updateBaktaUI();
  }
});

//gtdbtk
toggler("tax_class", "gtdbtk_dbDiv")

//checkm
toggler("run_checkm", "checkm_dbDiv")

//plasmids
toggler("run_plasmids", "genomad_dbDiv")