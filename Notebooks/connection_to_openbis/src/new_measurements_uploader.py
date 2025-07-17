from watchdog.observers import Observer
from watchdog.events import FileSystemEventHandler
import time
import os
import argparse
from pybis import Openbis

class NewFileHandler(FileSystemEventHandler):
    def __init__(self):
        self.measurements_folders = {}

    def on_created(self, event):
        if not event.is_directory:
            print(f"New file detected: {event.src_path}")
            self.custom_function(event.src_path)