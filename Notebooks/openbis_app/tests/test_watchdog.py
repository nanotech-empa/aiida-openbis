import time
from watchdog.observers import Observer
from watchdog.events import FileSystemEventHandler
import os
import utils

# Define the folder to watch
WATCHED_FOLDER = "/home/jovyan/aiida-openbis/Notebooks/connection_to_openbis/test_dataset_uploader"

# Create openBIS session
CONFIG = utils.read_json("config.json")
CONFIG_ELN = utils.get_aiidalab_eln_config()
DATA_MODEL = utils.read_yaml("/home/jovyan/aiida-openbis/Notebooks/Metadata_Schemas_LinkML/materialMLinfo.yaml")
OPENBIS_SESSION, SESSION_DATA = utils.connect_openbis(CONFIG_ELN["url"], CONFIG_ELN["token"])

class UploadHandler(FileSystemEventHandler):
    def on_created(self, event):
        if not event.is_directory:
            filename = os.path.basename(event.src_path)
            OPENBIS_SESSION.new_dataset(type = "RAW_DATA", sample = "20250602072805092-10288", files = [event.src_path]).save()
            print(f"📁 File uploaded: {filename}")

if __name__ == "__main__":
    print(f"🔍 Watching folder: {WATCHED_FOLDER}")
    event_handler = UploadHandler()
    observer = Observer()
    observer.schedule(event_handler, path="/home/jovyan/aiida-openbis/Notebooks/connection_to_openbis/test_dataset_uploader", recursive=True)
    observer.schedule(event_handler, path="/home/jovyan/aiida-openbis/Notebooks/connection_to_openbis/test_dataset_uploader_2", recursive=True)
    observer.start()

    try:
        while True:
            time.sleep(1)
    except KeyboardInterrupt:
        observer.stop()

    observer.join()