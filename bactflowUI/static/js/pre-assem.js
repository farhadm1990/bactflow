

const runButton = document.getElementById("run-bt");
const stopButton = document.getElementById("stop-bt");
const helpButton = document.getElementById("help-bt");
let eventSource = null;


const filePicker = document.getElementById("filePicker");
const fastqDirInput = document.getElementById("fastq_dir");


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
  const formData = new FormData(form);
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

// initiate concatenation
function catFastq(){
  const lsfq = document.getElementById("fastq-ls");
  const lsDiv = document.getElementById("fastq-ls-div");
  const form = document.getElementById('preForm');
  const spinFastqLs = document.getElementById("spin-fastq-list");
  const formData = new FormData(form);
  spinFastqLs.style.display = "block";
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
    lsfq.innerHTML = `<p>Error: ${error.message}</p>`;
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

// run read stat
function read_stat(){
  const spinSeq = document.getElementById("spin-seqkit");
  spinSeq.style.display = "block";
  const concatDiv = document.getElementById("concat_reads");
  form = document.getElementById('preForm');
  formData = new FormData(form);
  fetch("/reads-stat", {
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
 
      const seqDiv = document.getElementById("seq-div");
      const outputDiv = document.getElementById("output-seqkit");
      const threshIn = document.getElementById("threhold-input");
      
      
      if (data.html_output){
        spinSeq.style.display = "none";
        seqDiv.style.display = "block";
        threshIn.style.display = "block";
        concatDiv.value = "false";
        seqDiv.innerHTML = data.html_output;

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
          
        } else {
            console.error("❌ 'num_seqs' column not found in table.");
        }
      } else {
        seqDiv.innerHTML = `<p>Error: ${data.error}</p>`;
      }
      // data tables 
      setTimeout(() => {
        if ($.fn.DataTable) {
          if ($.fn.DataTable.isDataTable("#stats-table")) {
            table.DataTable().destroy();
        }
          $('#stats-table').DataTable({
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

    })
    .catch(error => {
      console.error("Error fetching read-stat:", error);
      document.getElementById("output-seqkit").innerHTML = `<p>Error: ${error.message}</p>`;
    });
}


// read visualization 
function read_vis(){
  const qualDiv = document.getElementById("qual-plot");
  const form = document.getElementById("preForm");
  const formData = new FormData(form);
  const spinVis = document.getElementById("spin-vis");
  qualDiv.style.display = "block";
  spinVis.style.display = "block";
  const qualHeatDiv = document.getElementById("qual-heat-plot");
  const qualHeatDivPool = document.getElementById("qual-heat-pool-plot");
  const range1Div = document.getElementById("range1Div");
  const range2Div = document.getElementById("range2Div");
  const range3Div = document.getElementById("range3Div");

  fetch("/plot-qual", {method: "POST", body: formData})
  .then(response => {
    if(!response.ok){
      return response.json().then(data => {
        throw new Error(data.error);
      });
    } else {
      return response.json()
    }
  })
  .then( data =>{
    if (data.error){
      qualDiv.innerHTML = `<p style="color: red;"><strong>Error:</strong> ${data.error}</p>`;
    }
    if(data.graph){
      spinVis.style.display = "none";
      qualHeatDiv.style.display = "block";
      qualHeatDivPool.style.display = "block";
      range1Div.style.display = "flex";
      range2Div.style.display = "flex";
      range3Div.style.display = "flex";

      let width1 = qualDiv.clientWidth;
      let width2 = qualHeatDiv.clientWidth;
      let width3 = qualHeatDivPool.clientWidth;
      let height1 = qualDiv.clientHeight;
      let height2 = qualHeatDiv.clientHeight;
      let height3 = qualHeatDivPool.clientHeight;

      document.getElementById("range1").addEventListener("input", function() {
        width1 = parseInt(this.value);
        document.getElementById("width1-value").textContent = width1;
        Plotly.relayout(qualDiv, {width: width1})
      });
      document.getElementById("range1H").addEventListener("input", function() {
        height1 = parseInt(this.value);
        document.getElementById("height1-value").textContent = height1;
        Plotly.relayout(qualDiv, {height: height1})
      });

      document.getElementById("range2").addEventListener("input", function() {
        width2 = parseInt(this.value);
        document.getElementById("width2-value").textContent = width2;
        Plotly.relayout(qualHeatDiv, {width: width2})
      });
      document.getElementById("range2H").addEventListener("input", function() {
        height2 = parseInt(this.value);
        document.getElementById("height2-value").textContent = height2;
        Plotly.relayout(qualHeatDiv, { height: height2})
      });

      document.getElementById("range3").addEventListener("input", function() {
        width3 = parseInt(this.value);
        document.getElementById("width3-value").textContent = width3;
        Plotly.relayout(qualHeatDivPool, {width: width3})
      });
      document.getElementById("range3H").addEventListener("input", function() {
        height3 = parseInt(this.value);
        document.getElementById("height3-value").textContent = height3;
        Plotly.relayout(qualHeatDivPool, { height: height3})
      });

      console.log("This is widht " + height2);
      var layout1 = {
        width: width1,   
        height: height1,                 
        // margin: { l: 5, r: 5, t: 5, b: 5 }
    };

      var layout2 = {
        width: width1,  
        height: height2,                 
        // margin: { l: 5, r: 5, t: 5, b: 5 }
    };

    var layout3 = {
      width: width3,  
      height: height3,                 
      // margin: { l: 0.5, r: 0.5, t: 0.5, b: 0.5 }
  };

    
    var plotData1 = JSON.parse(data.graph);
    var plotData2 = JSON.parse(data.graph_heat);
    var plotData3 = JSON.parse(data.graph_heat2);

    plotData1.layout = { ...plotData1.layout, ...layout1 };
    plotData2.layout = { ...plotData2.layout, ...layout2 };
    plotData3.layout = { ...plotData3.layout, ...layout3 };

    
    Plotly.newPlot(qualDiv, plotData1.data, plotData1.layout);
    Plotly.newPlot(qualHeatDiv, plotData2.data, plotData2.layout);
    Plotly.newPlot(qualHeatDivPool, plotData3.data, plotData3.layout);
    }
  })
  .catch(error => qualDiv.innerHTML = console.error("Error fetching plot", error));
}

const statBtn = document.getElementById('read-stat');
statBtn.addEventListener("click", read_stat
);

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

  document.getElementById("preForm").addEventListener("input", () =>{
    saveFormData();
  });
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

  toggleDiv();

  motherDiv.addEventListener("change", toggleDiv)
});
};


