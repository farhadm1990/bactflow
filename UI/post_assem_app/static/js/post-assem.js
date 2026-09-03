const outputDiv = document.getElementById("output-bactflow");
const runButton = document.getElementById("run-bt");
const stopButton = document.getElementById("stop-bt");
const helpButton = document.getElementById("help-bt");
let eventSource = null;
let streamGeneration = 0;
const STOP_USER_MSG = "Bactflow has been stopped by the user!";
// const form = document.getElementById("postForm");
// const formData = new FormData(form);
const filePicker = document.getElementById("filePicker");
const fastqDirInput = document.getElementById("fastq_dir");






// connect to stream when on assembly
function connectToStream(action){
  if (eventSource){
    console.log("Stream already connected");
    return;
  }

  window.__bactflowUserStopped = false;
  BactflowTerminal.init();
  const gen = ++streamGeneration;
  let logBuffer = [];
  let flushTimer = null;

  const flushLog = () => {
    flushTimer = null;
    if (gen !== streamGeneration || window.__bactflowUserStopped) {
      logBuffer = [];
      return;
    }
    if (!logBuffer.length) {
      return;
    }
    const lines = logBuffer.filter((line) => {
      const t = String(line || "").trim();
      return !/exited with code\s+-?\d+|Process completed|Stream disconnected/i.test(t);
    });
    logBuffer = [];
    if (lines.length) {
      BactflowTerminal.appendMany(lines);
    }
  };

  eventSource = new EventSource(`/stream_bactflow?action-assem=${action}`);
  if (window.BactflowMeters) {
    BactflowMeters.attachStream(eventSource);
  }
  
  eventSource.onmessage = (event) =>{
    if (gen !== streamGeneration || window.__bactflowUserStopped) {
      return;
    }
    logBuffer.push(event.data);
    if (!flushTimer) {
      flushTimer = setTimeout(flushLog, 200);
    }
  };

  eventSource.onerror = (error) =>{
    console.error("Error in streaming output:", error);
    if (flushTimer) {
      clearTimeout(flushTimer);
      flushTimer = null;
    }
    logBuffer = [];
    try {
      eventSource.onmessage = null;
      eventSource.onerror = null;
      eventSource.close();
    } catch (e) { /* ignore */ }
    if (eventSource) {
      eventSource = null;
    }
    if (gen !== streamGeneration || window.__bactflowUserStopped) {
      return;
    }
    BactflowProcessEta.finalizeAll(BactflowTerminal, "stopped");
    BactflowTerminal.append("Stream disconnected.", true);
    BactflowTerminal.setStatus("Stream disconnected", "warn");
    updateButtonStates("stopped");
    refreshCircularPlotButton();
  };
}

//disconnect function
function disconnectStream(){
  streamGeneration += 1;
  if(eventSource){
    try {
      eventSource.onmessage = null;
      eventSource.onerror = null;
    } catch (e) { /* ignore */ }
    eventSource.close();
    eventSource = null;
  }
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
  // switch actions
  switch (action){
    case "run":
      {
        
        fetch(`/run_bactflow?action-assem=${action}`, 
          { method: "POST", 
            body : formData }
          )
        .then(async (response) => {
          const message = await response.text();
          BactflowTerminal.clear();
          if (!response.ok) {
            BactflowTerminal.append(message || "Error starting BactFlow.", true);
            BactflowTerminal.setStatus("Failed to start", "error");
            document.getElementById('run-bt').disabled = false;
            document.getElementById('stop-bt').disabled = true;
            document.getElementById('help-bt').disabled = false;
            return;
          }

          BactflowTerminal.append("Bactflow started :)", true);
          BactflowTerminal.setStatus("Running...", "run");
          updateButtonStates("running");
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

      fetch(`/run_bactflow?action-assem=${action}`, { method: "POST" })
      .then((response) => {
        if(!response.ok) {
          BactflowTerminal.append("Error showing help for BactFlow. It might already be running?!", true);
          document.getElementById('run-bt').disabled = false;
          document.getElementById('help-bt').disabled = false;
          return;
        };
        BactflowTerminal.clear();
        BactflowTerminal.append("Bactflow's help menu", true);

        // const action = "help";
        updateButtonStates("running");
        connectToStream(action);
        
        
      })
      .catch((error) => {
        BactflowTerminal.append("Failed to give you BactFlow help!" + error.message, true);
        runButton.disabled = false;
        
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

    const action = (e.submitter && e.submitter.value) || new FormData(postForm).get("action-assem");
    if (!action) {
      return;
    }
    if (action !== "stop") {
      window.__bactflowUserStopped = false;
      BactflowTerminal.clear();
      document.getElementById("run-bt").disabled = true;
      document.getElementById("stop-bt").disabled = false;
      document.getElementById("help-bt").disabled = true;
    }

    run_wf(action);
  });

  // Direct Stop click — same strategy as assembly, guaranteed even if submitter is missing.
  const stopBt = document.getElementById("stop-bt");
  if (stopBt) {
    stopBt.addEventListener("click", (e) => {
      e.preventDefault();
      e.stopPropagation();
      run_wf("stop");
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
  if (childId === "genetype") {
    childDiv.style.display = visible ? "flex" : "none";
  } else {
    childDiv.style.display = visible ? "block" : "none";
  }
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
    const baktaSel = document.getElementById("bakta_annot");
    const circBtn = document.getElementById("circ-plot-bt");
    if (!baktaSel || baktaSel.value !== "true") {
      stopBaktaReadyPoll();
      return;
    }
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
  const circularDiv = document.getElementById("circular-div");
  const circBtn = document.getElementById("circ-plot-bt");
  const circHint = document.getElementById("circ-plot-hint");
  if (!baktaSel) {
    return;
  }

  const enabled = baktaSel.value === "true";
  showChildPanel("bakta_dbDiv", enabled);
  showChildPanel("genetype", enabled);
  if (circularDiv) {
    circularDiv.style.display = enabled ? "block" : "none";
  }

  if (!enabled) {
    stopBaktaReadyPoll();
    if (circBtn) {
      circBtn.disabled = true;
    }
    if (circHint) {
      circHint.textContent = "Enable Bakta annotation to use the circular plot.";
    }
    const circDiv = document.getElementById("circDiv");
    if (circDiv) {
      circDiv.style.display = "none";
    }
    updateStrainFinderFromFields();
    return;
  }

  if (circBtn) {
    circBtn.disabled = true;
  }
  if (circHint) {
    circHint.textContent = "Checking for .gbk/.gbff/.gff/.gff3 annotation files…";
  }
  updateStrainFinderFromFields();
  refreshCircularPlotButton();
  startBaktaReadyPoll();
}

async function refreshCircularPlotButton() {
  const baktaSel = document.getElementById("bakta_annot");
  const circBtn = document.getElementById("circ-plot-bt");
  const circHint = document.getElementById("circ-plot-hint");
  if (!baktaSel || baktaSel.value !== "true") {
    return;
  }

  if (circHint) {
    circHint.textContent = "Checking for annotation files…";
  }

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

// function to show quast

function showReport(){
  const form = document.getElementById("postForm");
  const formData = new FormData(form);
  Promise.all([
    fetch("/check-quast", {method: "POST", body: formData}).then(res => res.json()),
    fetch("/check-bakta", {method: "POST", body: formData}).then(res => res.json()),
    fetch("/circular", {method: "POST", body: formData}).then(res => res.json()),
    fetch("/taxa-report", {method: "POST", body: formData}).then(res => res.json())
  ])
  
    // quast and bakta
    // .then(([quastData, baktaData, circPlt, taxData]) => {
    .then(([quastData, baktaData, circPlt, taxData]) => {
      let quastDiv = document.getElementById("quastDiv");
      let baktaDiv = document.getElementById("baktaDiv");
      let circDiv  = document.getElementById("circDiv");
      let circSpin = document.getElementById("spin-circ");
      let taxDiv = document.getElementById("taxa_class");
  
      if(quastData.exists){
        console.log("this is quast data" + quastData);
        console.log("✅ Quast report found! Loading...");
        
        

        quastDiv.style.display = "block";
        quastReport();
        clearInterval(reportCheckInterval);
      } else {
        console.log("⏳ Waiting for Quast report...");
        
      
        quastDiv.style.display = "none";
      }

      // bakta
      if (baktaData.exists) {
        console.log("✅ Bakta report found! Loading...");
        baktaDiv.style.display = "block";
        clearInterval(reportCheckInterval);
        baktaReport();  
      } else {
        console.warn("⏳ Bakta report not available yet...");
        baktaDiv.style.display = "none";
      }

      refreshCircularPlotButton();

      // circular genome (display only if already generated)
      if (circPlt.plot) {
        circDiv.style.display = "block";
        circSpin.style.display = "none";
        showCircularPlot(circPlt.plot);
      }

      // taxonomy
      if(taxData.exists){
       
        taxDiv.style.display = "block";
        
        taxReport();
        clearInterval(reportCheckInterval);
      } else {
        console.warn("⏳ Taxonomy table not available yet...");
        taxDiv.style.display = "none";
      }
    })
    .catch(error => console.error("Error checking Quast report:", error));
}

async function taxReport() {
  const form = document.getElementById("postForm");
  const formData = new FormData(form);

  try {
    
    let taxRes = await fetch("/taxa-report", { method: "POST", body: formData });
    let taxDiv = document.getElementById("taxa_class");
   
    if (!taxRes.ok) {
      console.error("❌ Taxa report fetch failed:", taxRes.status, taxRes.statusText);
      throw new Error(`HTTP error! Status: ${taxRes.status}`);
    }

    let taxData = await taxRes.json();

    let taxautputDiv = document.getElementById("output-taxa");

   
    if (taxData.exists && taxData.taxa_table) {
      taxDiv.style.display = "block";
      taxautputDiv.innerHTML = taxData.taxa_table;
      
      setTimeout(() => {
        if ($.fn.DataTable) {
          console.log("✅ Initializing DataTable for taxa...");
          if ($.fn.DataTable.isDataTable("#taxa-tab")) {
            table.DataTable().destroy();
        }
          $("#taxa-tab").DataTable({
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
      console.warn("⚠️ No Taxa table found!");
      taxautputDiv.innerHTML = `<p>No taxa table data available.</p>`;
    }

  } catch (error) {
    console.error("❌ Error fetching taxa report:", error);
    alert("An error occurred while fetching the taxa report.");
  }
}

async function quastReport() {
  const form = document.getElementById("postForm");
  const formData = new FormData(form);
  
  try {
   
    let quastResponse = await fetch("/quast-report", { method: "POST", body: formData });
    if (!quastResponse.ok) throw new Error("Failed to fetch QUAST report");
    let quastBlob = await quastResponse.blob();
    let quastUrl = URL.createObjectURL(quastBlob);
    
    let quastOutputDiv = document.getElementById("output-quast");
    quastOutputDiv.innerHTML = ""; 
    quastOutputDiv.innerHTML = `<iframe src="${quastUrl}" style="width: 100%; height: 100%; border: none;"></iframe>`;
   
    // for contig
    let contigResponse = await fetch("/contig-report", { method: "POST", body: formData });
    if (!contigResponse.ok) throw new Error("Failed to fetch Contig report");
    let contigBlob = await contigResponse.blob();
    let contigUrl = URL.createObjectURL(contigBlob);
    
    let contigOutputDiv = document.getElementById("contig-quast");
    contigOutputDiv.innerHTML = ""; 
    contigOutputDiv.innerHTML = `<iframe src="${contigUrl}" style="width: 100%; height: 100%; border: none;"></iframe>`;
    
  } catch (error) {
    console.error("Error:", error);
    alert("An error occurred while fetching the reports.");
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
  const el = document.getElementById(id);
  if (el) {
    el.scrollIntoView({ behavior: "smooth", block: "start" });
  }
}

async function readJsonError(res, fallback) {
  try {
    const data = await res.json();
    return { data, error: data.error || data.message || fallback };
  } catch (e) {
    return { data: {}, error: fallback };
  }
}

function downloadCircularPlot() {
  const img = document.getElementById("circularImg");
  if (!img || !img.src) {
    return;
  }
  const link = document.createElement("a");
  link.href = img.src;
  link.download = "circular_plot.png";
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

function showCircularPlot(src) {
  const img = document.getElementById("circularImg");
  const bar = document.getElementById("circ-plot-bar");
  const typeLabel = document.getElementById("circ-type-label");
  if (typeLabel) {
    typeLabel.textContent = `(${selectedGeneTypes().join(", ")})`;
  }
  if (img && src) {
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
      circBtn.textContent = "Create circular plot";
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

var reportCheckInterval = setInterval(showReport, 5000);


// gene annotation 
const genetype = document.getElementById("geneType");
const vars = genetype.value;
genetype.addEventListener("change", ()=> {
  let selectedOptions = Array.from(genetype.selectedOptions).map(option => option.value)
  selectedOptions = selectedOptions.join(",")

});

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
  snpsDiv.scrollIntoView({ behavior: "smooth", block: "start" });

  fetch("/snp-finder", {method: "POST", body: formData})
  .then(async (res) => {
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
  })
  .then(data => {
    if (spin) {
      spin.style.display = "none";
    }
    if (data.exists) {
      hideSnpError();
      if (snpsOUT) {
        snpsOUT.innerHTML = "<p>✅ Your SNP files have been created successfully </p>";
      }
    } else {
      showSnpError(data.error || "Your SNP files were not created. Check the reference genome path and try again.");
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
  svsDiv.scrollIntoView({ behavior: "smooth", block: "start" });

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
  .then(data => {
    if (spin) {
      spin.style.display = "none";
    }
    if (data.exists) {
      hideSvsError();
      if (svsOUT) {
        svsOUT.innerHTML = "<p>Your variant calling files have been created successfully.</p>";
      }
    } else {
      showSvsError(data.error || "Variant calling files were not created. Check the reference genome path and try again.");
    }
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
    } else {
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
    } else {
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
