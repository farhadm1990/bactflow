const outputDiv = document.getElementById("output-bactflow");
const runButton = document.getElementById("run-bt");
const stopButton = document.getElementById("stop-bt");
const helpButton = document.getElementById("help-bt");
let eventSource = null;
// const form = document.getElementById("postForm");
// const formData = new FormData(form);
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
  const form = document.getElementById("postForm");
  const formData = new FormData(form);

  
  const setupOnly = document.getElementById("setup_only").value;
  
  

  if(clickedButtonValue === "run" && setupOnly === "false"){
    const percentage = 0;
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
          progDiv.innerHTML = `<strong>Analysis was successfully done!</strong>`;
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
          // progDiv.style.display = "block";
          progressBar.style.width = "0%";
          progressBar.setAttribute("aria-valuenow", "0");
          progressBar.innerText = "0%";
          // updateProgress();
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

document.getElementById("postForm").addEventListener("submit", (e) => {
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
    let form = document.getElementById("postForm");
    
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
  let form = document.getElementById("postForm");
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

  document.getElementById("postForm").addEventListener("input", () =>{
    saveFormData();
  });
});

// function to show quast
function showReport(){
  const form = document.getElementById("postForm");
  const formData = new FormData(form);
  Promise.all([
    fetch("/check-quast", {method: "POST", body: formData}).then(res => res.json()),
    fetch("/check-prok", {method: "POST", body: formData}).then(res => res.json()),
    fetch("/circular", {method: "POST", body: formData}).then(res => res.json()),
    fetch("/taxa-report", {method: "POST", body: formData}).then(res => res.json())
  ])
  
    // quast and proka
    .then(([quastData, prokData, circPlt, taxData]) => {
      let quastDiv = document.getElementById("quastDiv");
      let prokkDiv = document.getElementById("prokkaDiv");
      let circDiv  = document.getElementById("circDiv");
      let circSpin = document.getElementById("spin-circ");
      let taxDiv = document.getElementById("taxa_class");
  
      if(quastData.exists){
        // console.log("✅ Quast report found! Loading...");
        
        

        quastDiv.style.display = "block";
        quastReport();
        clearInterval(quastCheckInterval);
      } else {
        console.log("⏳ Waiting for Quast report...");
        
      
        quastDiv.style.display = "none";
      }

      // Prokka
      if (prokData.exists) {
        console.log("✅ Prokka report found! Loading...");
        prokkDiv.style.display = "block";
        clearInterval(quastCheckInterval);
        prokkaReport();  
      } else {
        console.warn("⏳ Prokka report not available yet...");
        prokkDiv.style.display = "none";
      }

      // circular genome
      if(circPlt.plot){
        circDiv.style.display = "block";
        circSpin.style.display = 'block';
        circReport();
      } else {
        console.warn("⏳ Circular plot not available yet...");
        circDiv.style.display = "none";
      }

      // taxonomy
      if(taxData.exists){
       
        taxDiv.style.display = "block";
        taxReport();
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

async function prokkaReport() {
  const form = document.getElementById("postForm");
  const formData = new FormData(form);

  try {
    
    let prokkaResponse = await fetch("/check-prok", { method: "POST", body: formData });

   
    if (!prokkaResponse.ok) {
      console.error("❌ Prokka report fetch failed:", prokkaResponse.status, prokkaResponse.statusText);
      throw new Error(`HTTP error! Status: ${prokkaResponse.status}`);
    }

    let prokkaData = await prokkaResponse.json();

    let prokkaOutputDiv = document.getElementById("output-prokk");

   
    if (prokkaData.exists && prokkaData.count_tab) {
  
      prokkaOutputDiv.innerHTML = prokkaData.count_tab;

      setTimeout(() => {
        if ($.fn.DataTable) {
          console.log("✅ Initializing DataTable...");
          if ($.fn.DataTable.isDataTable("#prokka-tab")) {
            table.DataTable().destroy();
        }
          $("#prokka-tab").DataTable({
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
      console.warn("⚠️ No Prokka table found!");
      prokkaOutputDiv.innerHTML = `<p>No gene count data available.</p>`;
    }

  } catch (error) {
    console.error("❌ Error fetching Prokka report:", error);
    alert("An error occurred while fetching the Prokka report.");
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
      console.error("❌ Prokka report fetch failed:", circRes.status, circRes.statusText);
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
    console.error("❌ Error fetching Prokka report:", error);
    alert("An error occurred while fetching the Prokka report.");
  }
}

let quastCheckInterval = setInterval(showReport, 5000);


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
        abundImgOut.src = abundData.plot_abund + "?t=" + new Date().getTime();
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
    
    // let prevImg = document.getElementById("prevImgDiv");
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
      
      
      
      
      
      
      setTimeout(() => {
        prevImgOut.src = prevData.plot_prev + "?t=" + new Date().getTime();
          
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

//bakta 
toggler("bakta_annot", childID = "bakta_dbDiv")

//genetype 
toggler("prok_annot", "genetype")
toggler("bakta_annot", "genetype")