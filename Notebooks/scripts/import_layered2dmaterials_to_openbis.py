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

# Get layered 2D materials from Excel file and upload them to openBIS

df = pd.read_excel("/home/jovyan/backup_aiida-openbis/Inventories/2Dmaterials_list.xlsx", header = 2)

available_layered2dmaterials_openbis = [layered2dmaterial.props.all()["$name"] for layered2dmaterial in session.get_samples(type = "LAYERED_2D_MATERIAL")]

for column, item in df.iterrows():
    
    layered2dmaterial_name = item["Acronym"]
    
    if layered2dmaterial_name not in available_layered2dmaterials_openbis:
        
        # For the moment I am going to skip the materials with no name. However, this should be solved in the future because every object must have a name.
        if is_nan(layered2dmaterial_name) == False:
            layered2dmaterial_metadata_dict = {"$name": layered2dmaterial_name}

            if is_nan(item["Name"]) == False:
                layered2dmaterial_metadata_dict["supplier_own_name"] = item["Name"]
            
            if is_nan(item["Comments"]) == False:
                layered2dmaterial_metadata_dict["comments"] = item["Comments"]
                
            if is_nan(item["Receiving Date"]) == False:
                layered2dmaterial_metadata_dict["receive_date"] = item["Receiving Date"].strftime('%Y-%m-%d')
            
            # Get all the information about the 2d layers
            layers_2d_material = {"material": item["Material"], "number_layers": item["Layers"]}
            
            if is_nan(item["Impurities"]) == False:
                layers_2d_material["dopants"] = item["Impurities"]
                
            layered2dmaterial_metadata_dict["layers_2d"] = json.dumps(layers_2d_material)
            
            # Get all the information about the substrates
            if is_nan(item["Substrate"]) == False:
                all_substrates = item["Substrate"].split("/")
                
                substrates_2d_material = []
                for substrate in all_substrates:
                    
                    if "(" in substrate:
                        pattern = r'(\w+)\s*\((.*?)\)' # Used to get the number of substrates which is inside parameters in the Excel file
                        matches = re.findall(pattern, substrate)

                        if matches:
                            substrates_2d_material.append({"material": matches[0][0], "thickness": matches[0][1]})
                        else:
                            substrates_2d_material.append({"material": substrate})
                
                if len(substrates_2d_material) > 0:
                    layered2dmaterial_metadata_dict["substrates"] = json.dumps(substrates_2d_material)
            
            if is_nan(item["Supplier"]) == False:
                if item["Supplier"].strip()!="":
                    layered2dmaterial_metadata_dict["supplier"] = {"$name": item["Supplier"]}
        
            layered2dmaterial_parents = []
                
            # Connect supplier with crystal if there is a supplier available
            if "supplier" in layered2dmaterial_metadata_dict:
                
                available_suppliers_openbis = [supplier.props.all()['$name'] for supplier in session.get_samples(type = "SUPPLIER")]
                
                # Create Supplier object or get it from openBIS (in case it is already there)
                if layered2dmaterial_metadata_dict["supplier"]["$name"] not in available_suppliers_openbis:
                    supplier_object = create_object_openbis(
                        session, 
                        "SUPPLIER", 
                        "/INSTITUTIONS/SUPPLIERS/SUPPLIER_COLLECTION", 
                        layered2dmaterial_metadata_dict["supplier"],
                        parents = []
                    )
                else:
                    supplier_object = get_openbis_object_by_property(session, 
                                                                    "SUPPLIER", 
                                                                    layered2dmaterial_metadata_dict["supplier"]["$name"]
                    )
                
                # Remove supplier information from the 2D layered material metadata because it should be a connection to a SUPPLIER object and not a parameter
                layered2dmaterial_metadata_dict.pop("supplier")
                
                # Append supplier parent
                layered2dmaterial_parents.append(supplier_object)
                    
            layered2dmaterial_object = create_object_openbis(
                session, 
                "2D_LAYERED_MATERIAL", 
                "/MATERIALS/2D_MATERIALS/2D_MATERIAL_COLLECTION", 
                layered2dmaterial_metadata_dict,
                parents = layered2dmaterial_parents
            )
            
        else:
            print(f"Layered 2D Material does not have a name")
            
    else:
        print(f"Layered 2D Material {layered2dmaterial_name} already in openBIS.")