import subprocess
filepath = ""
docker_container = "openbis-app-a"
docker_container_folder = "/data/openbis/eln-lims-dropbox/O+DEFAULT+DEFAULT+OBS17938+PUBLICATION_DATA+paper/"
subprocess.run(["docker", "exec", "openbis-app-a", "mkdir", docker_container_folder])
subprocess.run(["docker", "exec", "openbis-app-a", "chmod", "go+w", docker_container_folder])
subprocess.run(["docker", "cp", filepath, f"{docker_container}:{docker_container_folder}"])