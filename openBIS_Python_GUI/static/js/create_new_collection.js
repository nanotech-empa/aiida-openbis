document.getElementById('return-menu').addEventListener('click', function() {
    window.location.href = '/';
});

document.getElementById('save-new-project').addEventListener("click", function() {
        var newProjectSpace = document.getElementById("space-selector");
        var newProjectName = document.getElementById("project-name");
        var newProjectDescription = document.getElementById("project-description");
        var data = {
            projectSpace: newProjectSpace.value,
            projectName: newProjectName.value,
            projectDescription: newProjectDescription.value
        }
        var xhr = new XMLHttpRequest();
        xhr.open("POST", "/create_new_project", true);
        xhr.setRequestHeader("Content-Type", "application/json;charset=UTF-8");

        xhr.onload = function() {
            if (xhr.readyState === 4 && xhr.status === 200) {
                var response = JSON.parse(xhr.responseText);
                alert(response.result);
            }
        };
        xhr.send(JSON.stringify({ data: data }));
});