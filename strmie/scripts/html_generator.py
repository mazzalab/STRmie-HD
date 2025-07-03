import pandas as pd
import numpy as np
import json


string_html="""

<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Excel Viewer & Editor</title>
    <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/css/bootstrap.min.css">
    <link rel="stylesheet" href="https://cdn.datatables.net/1.11.5/css/jquery.dataTables.min.css">
    <link rel="stylesheet" href="https://cdn.datatables.net/colreorder/1.5.0/css/colReorder.dataTables.min.css">
    <script src="https://code.jquery.com/jquery-3.6.0.min.js"></script>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/xlsx/0.18.5/xlsx.full.min.js"></script>
    <script src="https://cdn.datatables.net/1.11.5/js/jquery.dataTables.min.js"></script>
    <script src="https://cdn.datatables.net/colreorder/1.5.0/js/dataTables.colReorder.min.js"></script>
    <script src="https://cdn.jsdelivr.net/npm/chart.js"></script>

    <style>
        .table-container {
            max-height: 250px;
            overflow: auto;
            white-space: nowrap;
        }

        thead {
            position: sticky;
            top: 0;
            background: white;
            z-index: 2;
        }

        #sampleDetails2 {
            margin-top: 20px;
        }

        #sampleDetails2 p {
            font-size: 16px;
        }

        #barChart {
            display: none;
            max-width: 90%;
            max-height: 800px;
            margin-top: 20px;
            border: 2px solid #ddd;
            box-shadow: 2px 2px 10px rgba(0,0,0,0.2);
        }

        #sampleDetails {
            margin-top: 20px;
        }

        #sampleDetails p {
            font-size: 16px;
        }

        #sampleSelect {
            max-height: 200px;
            overflow-y: auto;
        }

        .spacer {
            height: 250px;
        }

        .selected-row {
            background-color: rgba(0, 123, 255, 0.2) !important; /* Sfondo azzurro chiaro per la riga */
        }

        .highlighted-cell {
            background-color: rgba(255, 193, 7, 0.6) !important; /* Sfondo giallo per il Sample Name */
            font-weight: bold;
        }




    </style>
</head>
<body>
    <div class="container mt-4">
        <h2 class="text-center mb-4">Upload the strmie report (Excel file)</h2>
        
        <div class="mb-3">
            <input type="file" id="fileInput" class="form-control" accept=".xlsx, .xls">
        </div>

        <h3 class="text-center">Select the raw_counts folder containing data for bar plot generation</h3>
        <input type="file" id="folderInput" webkitdirectory class="form-control mt-3" />
        
        <div class="table-container border rounded p-2 mt-3">
            <table id="excelTable" class="display table table-striped table-bordered" style="width:100%">
                <thead></thead>
                <tbody></tbody>
            </table>
        </div>

        <div style="height: 40px;"></div>
        
        <!-- Dropdown to choose between CAG_repeats and CCG_repeats --> 
        <div class="mb-3">
            <label for="chartSelect"><strong>Select the chart to display, then reselect the row of the desired sample:</strong></label>
            <select id="chartSelect" class="form-control" style="position: relative;">
                <option value="CAG_repeats">CAG_repeats</option>
                <option value="CCG_repeats">CCG_repeats</option>
            </select>
        </div>
        
        <div id="sampleDetails2"></div>

        <div class="text-center">
            <canvas id="barChart"></canvas>
        </div>

<<<<<<< Updated upstream
        <button id="alertButtonCAA" type="button" class="btn btn-danger mt-3" style="display:none;">LOI CAA > 10</button>
        <button id="alertButtonCCA" type="button" class="btn btn-outline-danger" style="display:none;">LOI CCA > 10</button>
=======
        <button id="alertButtonCAA" type="button" class="btn btn-danger mt-3" style="display:none;">LOI CAA > 10%</button>
        <button id="alertButtonCCA" type="button" class="btn btn-outline-danger" style="display:none;">LOI CCA > 10%</button>
        <button id="alertButtonDOI" type="button" class="btn btn-warning mt-3" style="display:none;">DOI > 10%</button>
>>>>>>> Stashed changes

        <div id="sampleDetails"></div>

        <h3 class="text-center mt-4 fs-5">Generate a table for recalculating indices for samples whose peak values for Allele 1 and Allele 2 do not match those in the histogram</h3>

        <div class="mb-3">
            <label for="sampleSelect">Select Sample:</label>
            <select id="sampleSelect" class="form-control" style="position: relative;">
                <option value="">Select a Sample</option>
            </select>
        </div>

        <!-- Buttons for adding rows and saving the Excel -->
        <button id="addRow" class="btn btn-primary mb-3">Add Row</button>
        <button id="saveExcel" class="btn btn-success mb-3">Save as Excel</button>

        <div class="table-container border rounded p-2">
            <table id="newExcelTable" class="display table table-striped table-bordered" style="width:100%">
                <thead>
                    <tr>
                        <th>Sample</th>
                        <th>CAG_Allele_1</th>
                        <th>CAG_Allele_2</th>
                    </tr>
                </thead>
                <tbody></tbody>
            </table>
        </div>

        <div class="spacer"></div>
    </div>

    <script>
        let selectedFolder = null;
        let barChart = null;
        let currentChartType = "CAG_repeats"; // Default chart

        // Gestione della selezione della cartella
        document.getElementById('folderInput').addEventListener('change', function(event) {
            selectedFolder = event.target.files;
        });

        document.getElementById('fileInput').addEventListener('change', function (event) {
            let file = event.target.files[0];
            let reader = new FileReader();

            reader.onload = function (e) {
                let data = new Uint8Array(e.target.result);
                let workbook = XLSX.read(data, { type: 'array' });
                let firstSheet = workbook.SheetNames[0];
                let sheetData = XLSX.utils.sheet_to_json(workbook.Sheets[firstSheet], { header: 1 });

                if (sheetData.length === 0) return;

                let tableHead = document.querySelector("#excelTable thead");
                let tableBody = document.querySelector("#excelTable tbody");

                tableHead.innerHTML = "";
                tableBody.innerHTML = "";

                let headerRow = document.createElement("tr");
                sheetData[0].forEach(header => {
                    let th = document.createElement("th");
                    th.textContent = header;
                    headerRow.appendChild(th);
                });
                tableHead.appendChild(headerRow);

                let sampleIndex = sheetData[0].indexOf("Sample");
                let allele1Index = sheetData[0].indexOf("CAG_repeatsPeak_Allele_1");
                let allele2Index = sheetData[0].indexOf("CAG_repeatsPeak_Allele_2");

                if (sampleIndex === -1 || allele1Index === -1 || allele2Index === -1) {
                    alert("Error: One or more required columns (Sample, CAG_repeatsPeak_Allele_1, CAG_repeatsPeak_Allele_2) are missing in the Excel file!");
                    return;
                }

                let sampleValues = new Set();
                sheetData.slice(1).forEach(row => {
                    sampleValues.add(row[sampleIndex]);
                });

                let sampleSelect = document.getElementById('sampleSelect');
                sampleSelect.innerHTML = '<option value="">Select a Sample</option>';
                sampleValues.forEach(value => {
                    let option = document.createElement("option");
                    option.value = value;
                    option.textContent = value;
                    sampleSelect.appendChild(option);
                });

                sheetData.slice(1).forEach(row => {
                    let tr = document.createElement("tr");
                    row.forEach(cell => {
                        let td = document.createElement("td");
                        td.textContent = cell;
                        tr.appendChild(td);
                    });
                    tableBody.appendChild(tr);
                });

                let table = $('#excelTable').DataTable({
                    "scrollY": "200px",  
                    "scrollX": true,  
                    "paging": false,   
                    "info": false,     
                    "searching": false,
                    "colReorder": true  
                });

                $('#excelTable tbody').on('click', 'tr', function () {
                    let table = $('#excelTable').DataTable();
                    let rowData = table.row(this).data();
                    if (rowData) {
                        let sampleName = rowData[sampleIndex]?.trim();
                        let allele1Value = rowData[allele1Index]?.trim();
                        let allele2Value = rowData[allele2Index]?.trim();


 
                        // Rimuove evidenziazione da tutte le righe e celle
                        $('#excelTable tbody tr').removeClass('selected-row');
                        $('#excelTable tbody tr td').removeClass('highlighted-cell');
                        // Aggiunge classe alla riga selezionata
                        $(this).addClass('selected-row');
                        // Evidenzia la cella corrispondente al Sample Name
                        $(this).find('td').eq(sampleIndex).addClass('highlighted-cell');


                        if (!selectedFolder) {
                            alert("Please select a folder first.");
                            return;
                        }
                        
                        let table2FilePath = [...selectedFolder].find(file => file.name.includes(sampleName) && file.name.endsWith('.csv'));
                        if (!table2FilePath) {
                            alert("Table2 file not found for selected sample.");
                            return;
                        }

                        let reader2 = new FileReader();
                        reader2.onload = function(e) {
                            let csvData = e.target.result;
                            let rows = csvData.split("\\n").map(row => row.split(","));
                            let headers = rows[0];
                            let dataRows = rows.slice(1);

                            let currentChartIndex = headers.indexOf(currentChartType);
                            if (currentChartIndex === -1) {
                                alert(`Error: Column "${currentChartType}" not found in the CSV file.`);
                                return;
                            }

                            let repeatsData = dataRows.map(row => row[currentChartIndex]).filter(value => value != null);
                            let counts = repeatsData.reduce((acc, value) => {
                                acc[value] = (acc[value] || 0) + 1;
                                return acc;
                            }, {});




                            let sortedEntries = Object.entries(counts)
                                .map(([key, value]) => [parseInt(key), value])  // Converti le chiavi in numeri
                                .sort((a, b) => a[0] - b[0]);  // Ordina in base al valore numerico della chiave

                            let sortedLabels = sortedEntries.map(entry => entry[0]);  // Estrai le etichette ordinate
                            let sortedValues = sortedEntries.map(entry => entry[1]);  // Estrai i valori corrispondenti



                            

                            // Creare o aggiornare il grafico esistente
                            if (barChart) {
                                barChart.data.labels = sortedLabels;
                                barChart.data.datasets[0].data = sortedValues;
                                barChart.update();
                            } else {
                                let ctx = document.getElementById('barChart').getContext('2d');
                                barChart = new Chart(ctx, {
                                    type: 'bar',
                                    data: {
                                        labels: sortedLabels,  // Usa le etichette ordinate
                                        datasets: [{
                                            label: `${currentChartType} Occurrences`,
                                            data: sortedValues,  // Usa i valori ordinati
                                            backgroundColor: 'rgba(54, 162, 235, 0.5)',
                                            borderColor: 'rgba(54, 162, 235, 1)',
                                            borderWidth: 1
                                        }]
                                    },
                                    options: {
                                        responsive: true,
                                        scales: {
                                            y: {
                                                beginAtZero: true
                                            }
                                        }
                                    }
                                });
                            }


                            let loiCAAIndex = headers.indexOf("LOI_CAA");
                            let loiCCAIndex = headers.indexOf("LOI_CCA");
                            
                            let loiCAAValue = loiCAAIndex !== -1 ? parseFloat(rowData[loiCAAIndex]) || 0 : 0;
                            let loiCCAValue = loiCCAIndex !== -1 ? parseFloat(rowData[loiCCAIndex]) || 0 : 0;
                            
                            document.getElementById('alertButtonCAA').style.display = loiCAAValue > 10 ? 'block' : 'none';
                            document.getElementById('alertButtonCCA').style.display = loiCCAValue > 10 ? 'block' : 'none';

<<<<<<< Updated upstream
=======
                            // Recupera le intestazioni direttamente dalla tabella HTML
                            let table = $('#excelTable').DataTable();
                            let headerCells = table.columns().header().toArray();
                            let headerNames = headerCells.map(cell => cell.innerText.trim());
                            // Trova l’indice della colonna "DOI"
                            let doiIndex = headerNames.indexOf("DOI");
                            // Prendi il valore dalla riga selezionata
                            let doiRaw = rowData[doiIndex];
                            let doiValue = doiIndex !== -1 && !isNaN(parseFloat(doiRaw)) ? parseFloat(doiRaw.toString().trim()) : 0;
                            // Mostra o nasconde il pulsante
                            document.getElementById('alertButtonDOI').style.display = doiValue > 10 ? 'block' : 'none';


>>>>>>> Stashed changes
                            // Aggiungere i dettagli sotto il grafico
                            let details = `
                                <p><strong>Sample:</strong> ${sampleName}</p>
                                <p><strong>CAG Allele 1:</strong> ${allele1Value}</p>
                                <p><strong>CAG Allele 2:</strong> ${allele2Value}</p>
                            `;
                            let details2 = `
                                <p><strong>Sample:</strong> ${sampleName}</p>
                            `;
                            document.getElementById('sampleDetails').innerHTML = details;
                            document.getElementById('sampleDetails2').innerHTML = details2;

                        };
                        reader2.readAsText(table2FilePath);
                    }
                });
            };
            reader.readAsArrayBuffer(file);
        });

        document.getElementById("chartSelect").addEventListener("change", function() {
            currentChartType = this.value;
            if (barChart) {
                barChart.destroy();  // Distruggiamo il grafico esistente prima di crearne uno nuovo
                barChart = null;
            }
        });

        // Funzione per aggiungere una riga alla nuova tabella
        document.getElementById('addRow').addEventListener('click', function() {
            let selectedSample = $('#sampleSelect').val();
            if (!selectedSample) {
                alert("Please select a Sample first.");
                return;
            }
            $('#newExcelTable tbody').append('<tr><td contenteditable="true">' + selectedSample + '</td><td contenteditable="true"></td><td contenteditable="true"></td></tr>');
        });

        // Funzione per salvare la tabella modificata come Excel
        document.getElementById('saveExcel').addEventListener('click', function() {
            let wb = XLSX.utils.table_to_book(document.getElementById("newExcelTable"));
            XLSX.writeFile(wb, "CAG_data_for_recalculating_indices.xlsx");
        });
    </script>
</body>
</html>

"""


def create_html(outdir,path_imgs):
    with open(outdir+"/report.html", 'w') as f:
        f.write(string_html)
