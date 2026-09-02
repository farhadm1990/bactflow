const outputDiv = document.getElementById("output-bactflow");
const runButton = document.getElementById("run-bt");
const stopButton = document.getElementById("stop-bt");
const helpButton = document.getElementById("help-bt");
let eventSource = null;
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

  BactflowTerminal.init();
  let logBuffer = [];
  let flushTimer = null;

  const flushLog = () => {
    flushTimer = null;
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
    logBuffer.push(event.data);
    if (!flushTimer) {
      flushTimer = setTimeout(flushLog, 200);
    }
  };

  eventSource.onerror = (error) =>{
    console.error("Error in streaming output:", error);
    flushLog();
    BactflowProcessEta.finalizeAll(BactflowTerminal, "stopped");
    BactflowTerminal.append("Stream disconnected.", true);
    BactflowTerminal.setStatus("Stream disconnected", "warn");
    eventSource.close();
    eventSource = null;
    updateButtonStates("stopped");
    refreshCircularPlotButton();
  };
}

//disconnect function
function disconnectStream(){
  if(eventSource){
    eventSource.close();
    eventSource = null;
  }
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
        .then((response) => {
          BactflowTerminal.clear();
        if(!response.ok) {
          BactflowTerminal.append("Error starting BactFlow. It might already be running?!", true);
          BactflowTerminal.setStatus("Failed to start", "error");
          document.getElementById('run-bt').disabled = false;
          document.getElementById('help-bt').disabled = false;
          return;
        };
        

        BactflowTerminal.append("Bactflow started :)", true);
        BactflowTerminal.setStatus("Running...", "run");

      
        
        
        //now we start streatming here
       
        updateButtonStates("running");
        connectToStream(action);
      
      
    })
    
    .catch((error) => {
      BactflowTerminal.append("Failed to start BactFlow" + error.message, true);
      BactflowTerminal.setStatus("Failed to start", "error");
      
    });
        break;
      }
    
    case "stop":
      {
        document.getElementById('help-bt').disabled = false;
        document.getElementById('run-bt').disabled = false;
    


    fetch(`/run_bactflow?action-assem=${action}`, { 
      method: "POST"
    })
      .then((response) => response.text())
      .then((message) => {
        BactflowTerminal.clear();
        BactflowTerminal.append(message, true);
        BactflowTerminal.setStatus("Stopped", "warn");
        
        // disconnect
        disconnectStream();
        updateButtonStates('stopped')
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
    BactflowTerminal.clear();

    const action = e.submitter.value;

    document.getElementById("run-bt").disabled = true;
    document.getElementById("stop-bt").disabled = false;
    document.getElementById("help-bt").disabled = true;

    run_wf(action);
  });

  restoreFormData(); //from local storage

  postForm.addEventListener("input", () => {
    saveFormData();
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
    return;
  }

  if (circBtn) {
    circBtn.disabled = true;
  }
  if (circHint) {
    circHint.textContent = "Run BactFlow with Bakta enabled. The button activates when .gbk files are ready.";
  }
  refreshCircularPlotButton();
  startBaktaReadyPoll();
}

async function refreshCircularPlotButton() {
  const baktaSel = document.getElementById("bakta_annot");
  const circBtn = document.getElementById("circ-plot-bt");
  const circHint = document.getElementById("circ-plot-hint");
  if (!circBtn || !baktaSel || baktaSel.value !== "true") {
    return;
  }

  circBtn.disabled = true;
  if (circHint) {
    circHint.textContent = "Checking for Bakta results…";
  }

  const form = document.getElementById("postForm");
  const formData = new FormData(form);

  try {
    const res = await fetch("/check-bakta-ready", { method: "POST", body: formData });
    const data = await res.json();
    const ready = data.ready === true;
    circBtn.disabled = !ready;
    if (circHint) {
      circHint.textContent = ready
        ? `Bakta results found (${data.gbk_count || 0} .gbk file(s)). You can create the plot.`
        : (data.message || "Run BactFlow with Bakta enabled first.");
    }
    if (ready) {
      stopBaktaReadyPoll();
    }
  } catch (error) {
    circBtn.disabled = true;
    if (circHint) {
      circHint.textContent = "Could not check Bakta results.";
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
        document.getElementById("circularImg").src = circPlt.plot;
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

  formData.append("generate", "true");

  const circDiv = document.getElementById("circDiv");
  const circSpin = document.getElementById("spin-circ");
  const circImg = document.getElementById("circularImg");

  circDiv.style.display = "block";
  circSpin.style.display = "flex";
  circImg.style.display = "none";
  if (circBtn) {
    circBtn.disabled = true;
    circBtn.textContent = "Creating plot…";
  }

  try {
    const circRes = await fetch("/circular", { method: "POST", body: formData });
    if (!circRes.ok) {
      throw new Error(`HTTP ${circRes.status}`);
    }

    const circData = await circRes.json();
    circSpin.style.display = "none";

    if (circData.plot) {
      circImg.src = circData.plot;
      circImg.style.display = "block";
    } else if (circData.reason === "no_bakta_dir") {
      alert("Bakta annotation output not found. Run BactFlow with Bakta enabled first.");
    } else {
      alert(circData.error || "Could not create the circular plot.");
    }
  } catch (error) {
    circSpin.style.display = "none";
    console.error("Error creating circular plot:", error);
    alert("An error occurred while creating the circular plot.");
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
      circImg.src = circData.plot;

    
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
function runSNP(){
  let form = document.getElementById("postForm");
  let formData = new FormData(form);
  let snpsDiv = document.getElementById("snpsDiv");
  let snpsOUT = document.getElementById("output-snp");

  fetch("/snp-finder", {method: "POST", body: formData})
  .then(res => res.json())
  .then(data => {


    if (data.exists) {
      snpsDiv.style.display = "block";
      snpsOUT.innerHTML = "<p>✅ Your SNP files have been created successfully </p>";
    } else {
      snpsDiv.style.display = "block";
      snpsOUT.innerHTML = "<p>⚠️ Your SNP files don't exist, perhaps an error occurred! </p>";
    }
  })
  .catch((error) => {
    console.error("❌ Error fetching SNPs report:", error);
    alert("An error occurred while fetching the SNPs report.");
  });
}

// SVS
function runVCF(){
  let form = document.getElementById("postForm");
  let formData = new FormData(form);
  let svsDiv = document.getElementById("svs");
  let svsOUT = document.getElementById("output-svs");

  fetch("/svs-finder", {method: "POST", body: formData})
  .then(res => res.json())
  .then(data => {


    if (data.exists) {
      svsDiv.style.display = "block";
      svsOUT.innerHTML = "<p>✅ Your SVS files have been created successfully </p>";
    } else {
      svsDiv.style.display = "block";
      svsOUT.innerHTML = "<p>⚠️ Your SVS files don't exist, perhaps an error occurred! </p>";
    }
  })
  .catch((error) => {
    console.error("❌ Error fetching SVS report:", error);
    alert("An error occurred while fetching the SVS report.");
  });
}

// abund plot 
async function abundRun(){

  const form = document.getElementById("postForm");
  const formData = new FormData(form);
  let abundSpin = document.getElementById("spin-abund");
  abundSpin.style.display = "block";
  try {

    let abundRes = await fetch("/abund-run", { method: "POST", body: formData });

   
    if (!abundRes.ok) {
      console.error("❌ Abundance report fetch failed:", abundRes.status, abundRes.statusText);
      throw new Error(`HTTP error! Status: ${abundRes.status}`);
    }

    let abundData = await abundRes.json();
    
    let abundImg = document.getElementById("abundImgDiv");
    let abundImgOut = document.getElementById("abundImg");
    let abundDiv = document.getElementById("output-abund");
    let abDIV = document.getElementById("abundDiv");
    abundImg.style.display = "block";
    if (abundData.exists  && abundData.abund_table) {
      
      abDIV.style.display = "block";
      abundDiv.innerHTML = abundData.abund_table
      setTimeout(() => {
        if ($.fn.DataTable) {
          console.log("✅ Initializing abund DataTable...");
          if ($.fn.DataTable.isDataTable("#abund-tab")) {
            table.DataTable().destroy();
        }
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
        } else {
          console.warn("⚠️ DataTables is not loaded.");
        }
      }, 500);

     

    
    } else {
      console.warn("⚠️ No abundance data was found!");
      abundDiv.innerHTML = `<p>No Abundance data available.</p>`;
    }

    if (abundData.exists && abundData.plot_abund) {
      
      abundSpin.style.display = "none";
      
      
      
      
      setTimeout(() => {
        abundImgOut.src = abundData.plot_abund;
          console.log("✅ Image src updated:", abundImgOut.src);
      }, 500);
     
    } else {
      console.warn("⚠️ No abundance data was found!");
      abundDiv.innerHTML = `<p>No Abundance data available.</p>`;
    }

  } catch (error) {
    console.error("❌ Error fetching Abundance report:", error);
    alert("An error occurred while fetching the Abundance report.");
  }
}

// prev plot 
async function prevRun(){

  const form = document.getElementById("postForm");
  const formData = new FormData(form);
 

  try {

    let prevRes = await fetch("/prev-run", { method: "POST", body: formData });

   
    if (!prevRes.ok) {
      console.error("❌ Abundance report fetch failed:", prevRes.status, prevRes.statusText);
      throw new Error(`HTTP error! Status: ${prevRes.status}`);
    }

    let prevData = await prevRes.json();
    
   

    let prevImg = document.getElementById("prevImgDiv");
    let prevImgOut = document.getElementById("prevImg");
    let prevDiv = document.getElementById("output-prev");
    let prevDIV = document.getElementById("prevDiv");
    prevDIV.style.display = "block";
    if (prevData.exists  && prevData.prev_table) {
      
     
      prevDiv.innerHTML = prevData.prev_table
      setTimeout(() => {
        if ($.fn.DataTable) {
          console.log("✅ Initializing prev DataTable...");
          if ($.fn.DataTable.isDataTable("#prev-tab")) {
            table.DataTable().destroy();
        }
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
        } else {
          console.warn("⚠️ DataTables is not loaded.");
        }
      }, 500);

     

    
    } else {
      console.warn("⚠️ No prevalance data was found!");
      prevDiv.innerHTML = `<p>No prevalance data available.</p>`;
    }

    if (prevData.exists && prevData.plot_prev) {
      
      
      
      prevImg.style.display = 'block';
      
      
      setTimeout(() => {
        prevImgOut.src = prevData.plot_prev;
          
      }, 500);
     
    } else {
      console.warn("⚠️ No prevalance data was found!");
      prevDiv.innerHTML = `<p>No Prevalance data available.</p>`;
    }

  } catch (error) {
    console.error("❌ Error fetching prevalance report:", error);
    alert("An error occurred while fetching the prevalance report.");
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
