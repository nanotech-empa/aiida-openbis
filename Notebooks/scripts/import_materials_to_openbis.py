# Import libraries
import json
import re
import pandas as pd
import numpy as np
import tempfile
from pybis import Openbis
import warnings
import argparse
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

def get_crystal_properties(metadata_dict, item):
    
    if is_nan(item["Elog-Plate"]) == False:
        metadata_dict["sample_plate"] = item["Elog-Plate"]
    
    # Get crystal material and crystal face
    match = re.match(r"([A-Za-z]+)\((\d+)\)", item["Material"])
    
    if match and item["Material"][-1] == ")": # Last char is a parenthesis
        metadata_dict["material"] = match.group(1)
        metadata_dict["face"] = match.group(2)
    else:
        metadata_dict["material"] = item["Material"]
    
    if is_nan(item["Size"]) == False:
        
        # Find all numbers inside the string
        numbers = re.findall(r'\d+', item["Size"])
        crystal_diameter = numbers[0] if len(numbers) > 0 else None
        crystal_height = numbers[1] if len(numbers) > 1 else None
        
        # json.dumps is necessary to convert the python dictionary to json object
        metadata_dict["diameter"] = json.dumps({"value": crystal_diameter, "unit": "mm"})
        metadata_dict["height"] = json.dumps({"value": crystal_height, "unit": "mm"})
    
    if is_nan(item["Comment"]) == False:
        metadata_dict["comments"] = item["Comment"]
    
    if is_nan(item["Manufacturer"]) == False:
        if item["Manufacturer"].strip()!="":
            metadata_dict["supplier"] = {"$name": item["Manufacturer"]}
            
    return metadata_dict

def get_2d_layer_properties(metadata_dict, item):
    
    if is_nan(item["Name"]) == False:
        metadata_dict["supplier_own_name"] = item["Name"]

    if is_nan(item["Comments"]) == False:
        metadata_dict["comments"] = item["Comments"]

    if is_nan(item["Receiving Date"]) == False:
        metadata_dict["receive_date"] = item["Receiving Date"].strftime('%Y-%m-%d')

    # Get all the information about the 2d layers
    layers_2d_material = {"material": item["Material"], "number_layers": item["Layers"]}

    if is_nan(item["Impurities"]) == False:
        layers_2d_material["dopants"] = item["Impurities"]
        
    metadata_dict["layers_2d"] = json.dumps(layers_2d_material)

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
            metadata_dict["substrates"] = json.dumps(substrates_2d_material)

    if is_nan(item["Supplier"]) == False:
        if item["Supplier"].strip()!="":
            metadata_dict["supplier"] = {"$name": item["Supplier"]}
    
    return metadata_dict

# Set up openBIS environment
if __name__ == "__main__":
    # Read command line arguments
    
    parser = argparse.ArgumentParser(description = 'Setup the openBIS database (create objects types, collections, etc.).')

    # Define the arguments with flags
    parser.add_argument('-o', '--openbis_url', type=str, help='OpenBIS URL', default = 'https://local.openbis.ch/openbis')
    parser.add_argument('-u', '--openbis_user', type=str, help='OpenBIS User', default = 'admin')
    parser.add_argument('-pw', '--openbis_pw', type=str, help='OpenBIS Password', default = '123456789')
    parser.add_argument('-type', '--material_type', type=str, help='Material type (Crystal, Molecule, 2D-Layer Material, ...)', default = None)
    parser.add_argument('-xls', '--materials_filepath', type=str, help='Excel file with inventory', default = None)

    args = parser.parse_args()
    
    openbis_url = args.openbis_url
    openbis_user = args.openbis_user
    openbis_pw = args.openbis_pw
    material_type = args.material_type
    materials_filepath = args.materials_filepath
    
    if material_type:
        material_type = material_type.lower()
        if material_type in ["crystal", "2d-layer material"]:
            if materials_filepath:
                # Connect to openBIS
                session = log_in(bisurl=openbis_url, bisuser=openbis_user, bispasswd=openbis_pw)

                # Get materials from Excel file and upload them to openBIS
                if material_type == "2d-layer material":
                    df = pd.read_excel(materials_filepath, header = 2)
                    material_openbis_object_type = "2D_LAYERED_MATERIAL"
                elif material_type == "crystal":
                    df = pd.read_excel(materials_filepath, sheet_name = "Drawer")
                    material_openbis_object_type = "CRYSTAL"
                    
                available_materials_openbis = [material.props.all()["$name"] for material in session.get_samples(type = material_openbis_object_type)]

                for column, item in df.iterrows():
                    
                    if material_type == "2d-layer material":
                        material_name = item["Acronym"]
                    else:
                        material_name = item["Elog Name"]
                    
                    if material_name not in available_materials_openbis:
                        
                        # For the moment I am going to skip the materials with no name. However, this should be solved in the future because every object must have a name.
                        if is_nan(material_name) == False:
                            
                            material_metadata_dict = {"$name": material_name}

                            # 2D-Layer Material Case
                            if material_type == "2d-layer material":
                                material_metadata_dict = get_2d_layer_properties(material_metadata_dict, item)
                            else:
                                material_metadata_dict = get_crystal_properties(material_metadata_dict, item)
                        
                            # Connect supplier with crystal if there is a supplier available
                            if "supplier" in material_metadata_dict:
                                
                                available_suppliers_openbis = [supplier.props.all()['$name'] for supplier in session.get_samples(type = "SUPPLIER")]
                                
                                # Create Supplier object or get it from openBIS (in case it is already there)
                                if material_metadata_dict["supplier"]["$name"] not in available_suppliers_openbis:
                                    supplier_object = create_object_openbis(
                                        session, 
                                        "SUPPLIER", 
                                        "/INSTITUTIONS/SUPPLIERS/SUPPLIER_COLLECTION", 
                                        material_metadata_dict["supplier"],
                                        parents = []
                                    )
                                else:
                                    supplier_object = get_openbis_object_by_property(session,
                                                                                    "SUPPLIER", 
                                                                                    material_metadata_dict["supplier"]["$name"]
                                    )
                                
                                # Remove supplier information from the 2D layered material metadata because it should be a connection to a SUPPLIER object and not a parameter
                                material_metadata_dict.pop("supplier")
                                
                                # Append supplier object to the property list
                                material_metadata_dict['has_supplier'] = supplier_object.permId
                                    
                            layered2dmaterial_object = create_object_openbis(
                                session, 
                                material_openbis_object_type, 
                                f"/MATERIALS/{material_openbis_object_type}S/{material_openbis_object_type}_COLLECTION", 
                                material_metadata_dict,
                                parents = []
                            )
                            
                        else:
                            print(f"{material_type} does not have a name")
                            
                    else:
                        print(f"{material_type} {material_name} already in openBIS.")
        else:
            print("Unrecognised material types. Allowed material types: Crystal, 2D-Layer Material, ...")
    else:
        print("Specify the material type. Ex.: Crystal, Molecule, 2D-Layer Material, ...")