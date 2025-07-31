import json
from pathlib import Path
import IPython.display
import ipywidgets as ipw
from pybis import Openbis
import ipyfilechooser
import yaml
import pandas as pd
import copy
import datetime
import os
from IPython.display import display
import re
from collections import deque

def full_listdir(path):
    return [f"{path}{os.sep}{filepath}" for filepath in os.listdir(path)]

def upload_datasets(openbis_session, method_object, support_files_widget, dataset_type):
    for filename in support_files_widget.value:
        file_info = support_files_widget.value[filename]
        write_file(file_info['content'], filename)
        openbis_session.new_dataset(type = dataset_type, sample = method_object, files = [filename]).save()
        os.remove(filename)

def is_valid_json(string):
    try:
        obj = json.loads(string)
        if isinstance(obj, dict):
            return True  # Accept JSON objects
        if isinstance(obj, list) and all(isinstance(item, dict) for item in obj):
            return True  # Accept lists of JSON objects
        return False  # Reject primitive lists and other types
    except (ValueError, TypeError):
        return False

def get_aiidalab_eln_config():
    eln_config = Path.home() / ".aiidalab" / "aiidalab-eln-config.json"
    eln_config.parent.mkdir(
        parents=True, exist_ok=True
    )  # making sure that the folder exists.

    with open(eln_config) as file:
        config = json.load(file)
        eln_url = config["default"]
        eln_config = config[eln_url]
        eln_token = eln_config["token"]
        return {"url": eln_url, "token": eln_token}

def read_json(filename):
    with open(filename, 'r') as file:
        return json.load(file)

def read_yaml(filepath: str):
    return yaml.safe_load(open(filepath, 'r'))

def write_json(json_content, filename):
    with open(filename, 'w') as file:
        json.dump(json_content, file, indent=4)

def read_file(filename):
    return open(filename, "rb").read()

def connect_openbis_aiida():
    try:
        eln_config = Path.home() / ".aiidalab" / "aiidalab-eln-config.json"
        eln_config.parent.mkdir(parents=True, exist_ok=True)  # making sure that the folder exists.
        config = read_json(eln_config)
        eln_url = "https://local.openbis.ch"
        eln_token = config[eln_url]["token"]
        session_data = {"url": eln_url, "token": eln_token}
        openbis_session = Openbis(eln_url, verify_certificates = False)
        openbis_session.set_token(eln_token)
    except ValueError:
        openbis_session = None
        session_data = {}
    
    return openbis_session, session_data

def connect_openbis(eln_url, eln_token):
    try:
        session_data = {"url": eln_url, "token": eln_token}
        openbis_session = Openbis(eln_url, verify_certificates = False)
        openbis_session.set_token(eln_token)
    except ValueError:
        print("Session is no longer valid. Please check if the token is still valid.")
        openbis_session = None
        session_data = {}
    
    return openbis_session, session_data

def write_file(file_content, filename):
    with open(filename, 'wb') as f:  # Save file content
        f.write(file_content)

def get_next_collection_code(openbis_session, collection_type):
    collections = openbis_session.get_experiments(type = collection_type)
    collection_number = 1
    if collections:
        collection_number = max(int(collection.code.rsplit('_')[-1]) for collection in collections) + 1
    return f"{collection_type}_{collection_number}"

def create_openbis_dataset(openbis_session, **kwargs):
    openbis_ds = openbis_session.new_dataset(**kwargs)
    openbis_ds.save()

def create_openbis_object(openbis_session, **kwargs):
    openbis_object = openbis_session.new_object(**kwargs)
    openbis_object.save()
    return openbis_object

def create_openbis_collection(openbis_session, **kwargs):
    collection_type = kwargs.get("type", "")
    collection_code = kwargs.get("code", "")
    if collection_code == "":
        collection_code = get_next_collection_code(openbis_session, collection_type)
        kwargs["code"] = collection_code
        
    collection = openbis_session.new_collection(**kwargs)
    collection.save()
    return collection

def get_openbis_collections(openbis_session, **kwargs):
    return openbis_session.get_collections(**kwargs)

def get_openbis_projects(openbis_session, **kwargs):
    return openbis_session.get_projects(**kwargs)

def get_openbis_objects(openbis_session, **kwargs):
    return openbis_session.get_objects(**kwargs)

def get_openbis_object(openbis_session, **kwargs):
    return openbis_session.get_object(**kwargs)

def get_openbis_collection(openbis_session, **kwargs):
    return openbis_session.get_collection(**kwargs)

def get_openbis_property_type(openbis_session, **kwargs):
    return openbis_session.get_property_type(**kwargs)

def get_current_datetime():
    return datetime.datetime.now()

def convert_datetime_to_string(dt):
    return dt.strftime('%Y%m%d%H%M%S')

def remove_digits_from_string(string):
    return ''.join([i for i in string if not i.isdigit()])

def is_nan(var):
    """Function to verify whether it is a NaN."""
    return var != var

def find_first_atomistic_model(openbis_session, openbis_object, openbis_type):
    if openbis_type == "ATOMISTIC_MODEL":
        first_object = openbis_object
        
    object_parents_identifiers = openbis_object.parents
    for parent_identifier in object_parents_identifiers:
        parent_object = openbis_session.get_object(parent_identifier)
        if parent_object.type == "GEOMETRY_OPTIMISATION" or parent_object.type == "ATOMISTIC_MODEL":
            first_object = find_first_atomistic_model(openbis_session, parent_object, parent_object.type)
        
    return first_object

def get_last_sample_experiment(openbis_session, sample_object, last_sample_preparation_object):
    # Automatically select the experiment where the last sample preparation task was saved
    last_sample_experiment = openbis_session.get_experiment(last_sample_preparation_object.attrs.experiment)
    
    # Check if the sample is being used in a set of measurements that were created in another experiment after the preparation steps
    sample_children = sample_object.get_children()
    for child in sample_children:
        if child.registrationDate > last_sample_preparation_object.registrationDate:
            last_sample_experiment = child.experiment
    
    return last_sample_experiment

def get_last_sample_instrument(openbis_session, last_sample_preparation_object):
    # Automatically select the instrument used in the last sample preparation task
    for parent in last_sample_preparation_object.parents:
        parent_object = openbis_session.get_object(parent)
        if parent_object.type == "INSTRUMENT":
            instrument_object = parent_object
            break
    
    return instrument_object

def is_numeric(s):
    try:
        float(s)  # Attempt conversion to float
        return True
    except ValueError:
        return False

def stringify_quantity_value(json_obj, unit_type):
    if isinstance(json_obj, str):
        json_obj = json.loads(json_obj)
        
    value = json_obj["value"]
    unit = json_obj[unit_type]
    quantity_value_str = ""
    if value:
        quantity_value_str = f"{value} {unit}"
        
    return quantity_value_str