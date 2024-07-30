# Import libraries
import json
import re
import pandas as pd
import numpy as np
import tempfile
from pybis import Openbis
import warnings
warnings.filterwarnings("ignore")

# Functions
def is_nan(value):
    return value != value

def read_file(filepath: str) -> str:
    with open(filepath, "rb") as f:
        return f.read()

def log_in(bisurl='openbis', bisuser='admin', bispasswd='changeit'):
    """Function to login to openBIS."""
    if Openbis(bisurl, verify_certificates=False).is_token_valid():
        session = Openbis(bisurl, verify_certificates=False)
    else:
        Openbis(bisurl, verify_certificates=False).login(bisuser, bispasswd, save_token=True)
        session = Openbis(bisurl, verify_certificates=False)
    return session

def create_object_openbis(session, object_type, object_experiment, object_props, parents):
    
    if len(parents) > 0:
        openbis_object = session.new_sample(
            type = object_type,
            experiment = object_experiment,
            props = object_props,
            parents = parents
        )
    else:
        openbis_object = session.new_sample(
            type = object_type,
            experiment = object_experiment,
            props = object_props
        )
    
    openbis_object.save()
    
    return openbis_object

def get_openbis_object_by_property(session, object_type, property_value):
    for object in session.get_samples(type = object_type):
        if object.props.all()['$name'] == property_value:
            return object

# Connect to openBIS
session = log_in(bisurl='openbis', bisuser='admin', bispasswd='123456789')

## Get crystals from Excel file and upload them to openBIS
df = pd.read_excel("/home/jovyan/backup_aiida-openbis/Inventories/Crystal.xlsx", sheet_name = "Drawer")

available_crystals_openbis = [crystal.props.all()["$name"] for crystal in session.get_samples(type = "CRYSTAL")]

for column, item in df.iterrows():
    
    crystal_name = item["Elog Name"]
    
    if crystal_name not in available_crystals_openbis:
        
        crystal_metadata_dict = {"$name": crystal_name}

        if is_nan(item["Elog-Plate"]) == False:
            crystal_metadata_dict["sample_plate"] = item["Elog-Plate"]
        
        # Get crystal material and crystal face
        match = re.match(r"([A-Za-z]+)\((\d+)\)", item["Material"])
        
        if match and item["Material"][-1] == ")": # Last char is a parenthesis
            crystal_metadata_dict["material"] = match.group(1)
            crystal_metadata_dict["face"] = match.group(2)
        else:
            crystal_metadata_dict["material"] = item["Material"]
        
        if is_nan(item["Size"]) == False:
            
            # Find all numbers inside the string
            numbers = re.findall(r'\d+', item["Size"])
            crystal_diameter = numbers[0] if len(numbers) > 0 else None
            crystal_height = numbers[1] if len(numbers) > 1 else None
            
            # json.dumps is necessary to convert the python dictionary to json object
            crystal_metadata_dict["diameter"] = json.dumps({"value": crystal_diameter, "unit": "mm"})
            crystal_metadata_dict["height"] = json.dumps({"value": crystal_height, "unit": "mm"})
        
        if is_nan(item["Comment"]) == False:
            crystal_metadata_dict["comments"] = item["Comment"]
        
        if is_nan(item["Manufacturer"]) == False:
            if item["Manufacturer"].strip()!="":
                crystal_metadata_dict["supplier"] = {"$name": item["Manufacturer"]}
        
        crystal_parents = []
            
        # Connect supplier with crystal if there is a supplier available
        if "supplier" in crystal_metadata_dict:
            
            available_suppliers_openbis = [supplier.props.all()['$name'] for supplier in session.get_samples(type = "SUPPLIER")]
            
            # Create Supplier object or get it from openBIS (in case it is already there)
            if crystal_metadata_dict["supplier"]["$name"] not in available_suppliers_openbis:
                supplier_object = create_object_openbis(
                    session, 
                    "SUPPLIER", 
                    "/INSTITUTIONS/SUPPLIERS/SUPPLIER_COLLECTION", 
                    crystal_metadata_dict["supplier"],
                    parents = []
                )
            else:
                supplier_object = get_openbis_object_by_property(session, 
                                                                "SUPPLIER", 
                                                                crystal_metadata_dict["supplier"]["$name"]
                )
            
            # Remove supplier information from the crystal metadata because it should be a connection to a SUPPLIER object and not a parameter
            crystal_metadata_dict.pop("supplier")
            
            # Append supplier parent
            crystal_parents.append(supplier_object)
                
        crystal_object = create_object_openbis(
            session, 
            "CRYSTAL", 
            "/MATERIALS/CRYSTALS/CRYSTAL_COLLECTION", 
            crystal_metadata_dict,
            parents = crystal_parents
        )
            
    else:
        print(f"Crystal {crystal_name} already in openBIS.")