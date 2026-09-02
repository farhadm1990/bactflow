const outputDiv = document.getElementById("output-bactflow");
const runButton = document.getElementById("run-bt");
const stopButton = document.getElementById("stop-bt");
const helpButton = document.getElementById("help-bt");
let eventSource = null;


const fastqDirInput = document.getElementById("fastq_dir");



//progress function
// getting values of the run buttons 
// let clickedButtonValue = null;
// document.addEventListener("click", function(event) {
//   if (event.target.name === "action-assem") {
//        clickedButtonValue = event.target.value;
     
//   }
// });

// const updateProgress = () =>{
//   const form = document.getElementById("runForm");
//   const formData = new FormData(form);

  
//   const setupOnly = document.getElementById("setup_only").value;
  
  

//   if(clickedButtonValue === "run" && setupOnly === "false"){
    
//     fetch("/progress", { method: "POST", body: formData } )
//     .then(response => response.json())
//     .then(data =>{
//       let percentage = data.completed;
//       let progressBar = document.querySelector(".progress-bar");
//       let progDiv = document.querySelector(".progress");
      
//       progDiv.style.display = "block";
      
//       progressBar.style.width = percentage + "%";
//       progressBar.setAttribute("aria-valuenow", percentage);

//       progressBar.innerText = Math.round(percentage) + "%";
      
//       if(percentage < 100){
//         setTimeout(updateProgress,2000);
//       } else {
//         progDiv.style.transition = "background-color 1s ease, color 1s ease, text-align 0s ease";
//         setTimeout(() => {
//           progDiv.style.backgroundColor = "green";
//           progDiv.style.color = "white";
//           progressBar.style.opacity = "0"; 
  
//           setTimeout(() => {
//               progressBar.style.display = "none"; 
//           }, 5);
  
        
//           progDiv.style.opacity = "5"; 
//           progDiv.innerHTML = `<strong>Assembly was successfully done!</strong>`;
//           progDiv.style.textAlign = "justify";
//           progDiv.style.textAlignLast = "center";
//           progDiv.style.height = '30px';
//           progDiv.style.padding = '5px';

  
//           setTimeout(() => {
//               progDiv.style.opacity = "1"; // Fade in the text
//           }, 20); // Small delay for smooth appearance
  
//       }, 200); // Slight delay after progress completes
//         // progDiv.style.backgroundColor = "green";
//         // progDiv.style.color = "white";
//         // progressBar.style.display = "none";
//         // progDiv.style.textAlign = "justify";
//         // progDiv.style.textAlignLast= "center"; 


//         progDiv.innerHTML = `<strong>Assembly was successfully done!</strong>`

//       }
//     })
//     .catch(error => console.error("Error fetching progress:", error))
//   } 
  
// }


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
  showEl(flyeCoverage, isFlye || isPacbio);
  showEl(nanofilterGroup, isLongRead);
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
  // switch actions
  switch (action ){
    case "run":
      {
        
        fetch(`/run_bactflow?action-assem=${action}`, { method: "POST", body : formData })
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

        // getting value of the setup only field
        let setOnly = document.getElementById("setup_only").value;
        // let progressBar = document.querySelector(".progress-bar");
        // let progDiv = document.querySelector(".progress");
          // progressBar.style.width = "0%";
          // progressBar.setAttribute("aria-valuenow", "0");
          // progressBar.innerText = "0%";
        // if(setOnly === "false"){
          // progDiv.style.display = "block";
          // progressBar.style.width = "0%";
          // progressBar.setAttribute("aria-valuenow", "0");
          // progressBar.innerText = "0%";
          // updateProgress();
        // } else {
        //   progDiv.style.display = "none";
        // };
        
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

document.getElementById("runForm").addEventListener("submit", (e) => {
 e.preventDefault();
 BactflowTerminal.clear();

  const action = e.submitter.value;
  const selectedAssembler = document.getElementById("assemblerDropdown").value;
  const longReads = document.getElementById("fastq_dir").value.trim();
  if (action === "run" && selectedAssembler === "unicycler") {
    const shortReads = document.getElementById("short_read_dir").value.trim();
    if (!longReads || !shortReads) {
      BactflowTerminal.append("Unicycler hybrid assembly needs both a long-read path and a short-read path.", true);
      document.getElementById('run-bt').disabled = false;
      document.getElementById('stop-bt').disabled = true;
      document.getElementById('help-bt').disabled = false;
      return;
    }
  }
  if (action === "run" && selectedAssembler === "spades" && !longReads) {
    BactflowTerminal.append("SPAdes needs an Illumina paired-end FASTQ directory (sample_R1 / sample_R2).", true);
    document.getElementById('run-bt').disabled = false;
    document.getElementById('stop-bt').disabled = true;
    document.getElementById('help-bt').disabled = false;
    return;
  }

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
    let form = document.getElementById("runForm");
    
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
  let form = document.getElementById("runForm");
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
  updateAssemblerUI();

  document.getElementById("runForm").addEventListener("input", () =>{
    saveFormData();
  });
});

// function to show quast
function checkForQuastReport(){
  const form = document.getElementById("runForm");
  const formData = new FormData(form);

  fetch("/check-quast", {method: "POST", body: formData})
    .then(response => response.json())
    .then(data => {
      let helpBtn = document.getElementById("help-bt");

      if(data.exists){
        console.log("✅ Quast report found! Loading...");
        
        let quastDiv = document.getElementById("quastDiv");
        quastDiv.style.display = "block";
        quastReport();
        clearInterval(quastCheckInterval);
      } else {
        console.log("⏳ Waiting for Quast report...");
        let quastDiv = document.getElementById("quastDiv");
      
        quastDiv.style.display = "none";
      }
    })
    .catch(error => console.error("Error checking Quast report:", error));
}


async function quastReport() {
  const form = document.getElementById("runForm");
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


let quastCheckInterval = setInterval(checkForQuastReport, 5000);




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
  fastqDirInput.addEventListener("change", (event) => {
    const selectedFiles = Array.from(event.target.files || []);

    if (selectedFiles.length > 0) {
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

//gtdbtk
toggler("tax_class", "gtdbtk_dbDiv")

//checkm
toggler("run_checkm", "checkm_dbDiv")