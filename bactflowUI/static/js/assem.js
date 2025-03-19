const outputDiv = document.getElementById("output-bactflow");
const runButton = document.getElementById("run-bt");
const stopButton = document.getElementById("stop-bt");
const helpButton = document.getElementById("help-bt");
let eventSource = null;


const filePicker = document.getElementById("filePicker");
const fastqDirInput = document.getElementById("fastq_dir");



//progress function
// getting values of the run buttons 
let clickedButtonValue = null;
document.addEventListener("click", function(event) {
  if (event.target.name === "action-assem") {
       clickedButtonValue = event.target.value;
     
  }
});

const updateProgress = () =>{
  const form = document.getElementById("runForm");
  const formData = new FormData(form);

  
  const setupOnly = document.getElementById("setup_only").value;
  
  

  if(clickedButtonValue === "run" && setupOnly === "false"){
    
    fetch("/progress", { method: "POST", body: formData } )
    .then(response => response.json())
    .then(data =>{
      let percentage = data.completed;
      let progressBar = document.querySelector(".progress-bar");
      let progDiv = document.querySelector(".progress");
      
      progDiv.style.display = "block";
      
      progressBar.style.width = percentage + "%";
      progressBar.setAttribute("aria-valuenow", percentage);

      progressBar.innerText = Math.round(percentage) + "%";
      
      if(percentage < 100){
        setTimeout(updateProgress,2000);
      } else {
        progDiv.style.transition = "background-color 1s ease, color 1s ease, text-align 0s ease";
        setTimeout(() => {
          progDiv.style.backgroundColor = "green";
          progDiv.style.color = "white";
          progressBar.style.opacity = "0"; 
  
          setTimeout(() => {
              progressBar.style.display = "none"; 
          }, 5);
  
        
          progDiv.style.opacity = "5"; 
          progDiv.innerHTML = `<strong>Assembly was successfully done!</strong>`;
          progDiv.style.textAlign = "justify";
          progDiv.style.textAlignLast = "center";
          progDiv.style.height = '30px';
          progDiv.style.padding = '5px';

  
          setTimeout(() => {
              progDiv.style.opacity = "1"; // Fade in the text
          }, 20); // Small delay for smooth appearance
  
      }, 200); // Slight delay after progress completes
        // progDiv.style.backgroundColor = "green";
        // progDiv.style.color = "white";
        // progressBar.style.display = "none";
        // progDiv.style.textAlign = "justify";
        // progDiv.style.textAlignLast= "center"; 


        progDiv.innerHTML = `<strong>Assembly was successfully done!</strong>`

      }
    })
    .catch(error => console.error("Error fetching progress:", error))
  } 
  
}


// connect to stream when on assembly
function connectToStream(action){
  if (eventSource){
    console.log("Stream already connected");
    return;
  }

  
  outputDiv.innerHTML = "";
  eventSource = new EventSource(`/stream_bactflow?action-assem=${action}`);
  
  eventSource.onmessage = (event) =>{
    
    outputDiv.innerHTML += event.data + '\n';
    outputDiv.scrollTop = outputDiv.scrollHeight;
    
  };

  eventSource.onerror = (error) =>{
    console.error("Error in streaming output:", error);
    outputDiv.innerHTML += "Stream disconnected.  \n";
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
  var selectedAssembler = this.value;

  var assemblers = ["flye", "unicycler", "megahit", "spades"];
  assemblers.forEach(assembler => {
    document.getElementById("run_" + assembler).value = "false";
  });
  if(selectedAssembler !=="none") {
    document.getElementById("run_" + selectedAssembler).value = "true";
  }
});

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
          outputDiv.innerHTML = "";
        if(!response.ok) {
          outputDiv.innerHTML += "Error starting BactFlow. It might already be running?!\n";
          document.getElementById('run-bt').disabled = false;
          document.getElementById('help-bt').disabled = false;
          return;
        };
        

        outputDiv.innerHTML += "Bactflow started :)\n";

        // getting value of the setup only field
        let setOnly = document.getElementById("setup_only").value;
        let progressBar = document.querySelector(".progress-bar");
        let progDiv = document.querySelector(".progress");
          progressBar.style.width = "0%";
          progressBar.setAttribute("aria-valuenow", "0");
          progressBar.innerText = "0%";
        if(setOnly === "false"){
          progDiv.style.display = "block";
          progressBar.style.width = "0%";
          progressBar.setAttribute("aria-valuenow", "0");
          progressBar.innerText = "0%";
          updateProgress();
        } else {
          progDiv.style.display = "none";
        };
        
        //now we start streatming here
       
        updateButtonStates("running");
        connectToStream(action);
      
      
    })
    
    .catch((error) => {
      outputDiv.innerHTML += "Failed to start BactFlow" + error.message + "\n";
      
    });
        break;
      }
    
    case "stop":
      {
        document.getElementById('help-bt').disabled = false;
        document.getElementById('run-bt').disabled = false;
        let progDiv = document.querySelector(".progress");
        progDiv.style.display = "none";


    fetch(`/run_bactflow?action-assem=${action}`, { 
      method: "POST"
    })
      .then((response) => response.text())
      .then((message) => {
        outputDiv.innerHTML = "";
        outputDiv.innerHTML += message + "\n";
        
        // disconnect
        disconnectStream();
        updateButtonStates('stopped')
      });
          break;
        }
        
      case "help":
      {
        outputDiv.innerHTML = "";
        document.getElementById('help-bt').disabled = true;
        let progDiv = document.querySelector(".progress");
        progDiv.style.display = "none";

        // hide quast
        const quastDiv = document.getElementById("quastDiv");
        quastDiv.style.display = "none";

      fetch(`/run_bactflow?action-assem=${action}`, { method: "POST" })
      .then((response) => {
        if(!response.ok) {
          outputDiv.innerHTML += "Error showing help for BactFlow. It might already be running?!\n";
          document.getElementById('run-bt').disabled = false;
          document.getElementById('help-bt').disabled = false;
          return;
        };
        outputDiv.innerHTML = "";
        outputDiv.innerHTML += "Bactflow's help menue!\n";

        // const action = "help";
        updateButtonStates("running");
        connectToStream(action);
        
        
      })
      .catch((error) => {
        outputDiv.innerHTML += "Failed to give you BactFlow help!" + error.message + "\n";
        runButton.disabled = false;
        
      });
        break;
      }
      
    
  }
}

document.getElementById("runForm").addEventListener("submit", (e) => {
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

//medaka basecaller model toggle
toggler("medaka_polish", "medakaModel");

//coverage parameter toggle
toggler("coverage_filter", "covDiv");

// filter read
toggler(motherId = "nanofilter", childId = "qualDiv");

// medaka tensor
toggler(motherId = "medaka_polish", childID = "tesnDiv")

//gtdbtk
toggler("tax_class", "gtdbtk_dbDiv")

//checkm
toggler("run_checkm", "checkm_dbDiv")