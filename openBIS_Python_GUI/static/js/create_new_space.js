document.getElementById('return-menu').addEventListener('click', function() {
    window.location.href = '/';
});

document.getElementById('save-new-space').addEventListener("click", function() {
        var newSpaceName = document.getElementById("space-name");
        var newSpaceDescription = document.getElementById("space-description");
        var data = {
            spaceName: newSpaceName.value,
            spaceDescription: newSpaceDescription.value
        }
        var xhr = new XMLHttpRequest();
        xhr.open("POST", "/create_new_space", true);
        xhr.setRequestHeader("Content-Type", "application/json;charset=UTF-8");

        xhr.onload = function() {
            if (xhr.readyState === 4 && xhr.status === 200) {
                var response = JSON.parse(xhr.responseText);
                alert(response.result);
            }
        };
        xhr.send(JSON.stringify({ data: data }));
});