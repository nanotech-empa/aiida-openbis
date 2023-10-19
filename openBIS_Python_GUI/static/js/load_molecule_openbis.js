function handleEnter(event) {
    if (event.key === "Enter") {
        // Call a function to send the entered data to Python
        sendDataToPython(document.getElementById("complex-input").value);
    }
}

function sendDataToPython(data) {
    // Use AJAX to send the data to your Python backend
    var xhr = new XMLHttpRequest();
    xhr.open("POST", "/molecule/update_molecule", true);
    xhr.setRequestHeader("Content-Type", "application/json;charset=UTF-8");
    xhr.onreadystatechange = function () {
        if (xhr.readyState === 4 && xhr.status === 200) {
            var response = JSON.parse(xhr.responseText);
            document.getElementById('molecule-id').value = response.moleculeId;
            document.getElementById('complex-input').value = response.moleculeEmpaNumber;
            document.getElementById('text-selector').value = response.moleculeBatch;
            document.getElementById('molecule-acronym').value = response.moleculeAcronym;
            document.getElementById('molecule-iupac').value = response.moleculeIupacName;
            document.getElementById('molecule-formula').value = response.moleculeFormula;
            document.getElementById('molecule-smiles').value = response.moleculeSmiles;
            document.getElementById('molecule-cas-number').value = response.moleculeCasNumber;
            document.getElementById('hazardous-checkbox').checked = response.moleculeHazardous;
            document.getElementById('molecule-hazardous-specify').value = response.moleculeHazardousSpecify;
            document.getElementById('light-checkbox').checked = response.moleculeLight;
            document.getElementById('oxygen-checkbox').checked = response.moleculeOxygen;
            document.getElementById('fridge-checkbox').checked = response.moleculeFreezer;
            document.getElementById('dry-checkbox').checked = response.moleculeDry;
            document.getElementById('other-condition-checkbox').checked = response.moleculeStorageConditionOther;
            document.getElementById('other-condition-text').value = response.moleculeStorageConditionOtherSpecify;
            document.getElementById('molecule-supplier-amount').value = response.moleculeAmount;
            document.getElementById('molecule-receiving-date').value = response.moleculeReceivingDate;
            document.getElementById('molecule-comments').value = response.moleculeComments;
        }
    };
    xhr.send(JSON.stringify({ data: data }));
}

document.getElementById('return-menu').addEventListener('click', function() {
    window.location.href = '/';
});