# Import libraries
import json
import re
import os
import rdkit
from rdkit.Chem import AllChem, Draw, rdMolDescriptors
import pandas as pd
import numpy as np
import tempfile
from pybis import Openbis
import warnings
import argparse
warnings.filterwarnings("ignore")

# Global variables
MATERIALS_TYPES_LIST = ["crystal", "2d-layer material", "molecule"]

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

def get_full_filepath(folderpath: str) -> list[str]:
    return [f"{folderpath}/{file}" for file in os.listdir(folderpath) if file.endswith(".cdxml")]

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

def create_dataset_openbis(session, dataset_type, object_id, dataset_filepath):
    raw_ds = session.new_dataset(type = dataset_type, object = object_id, file = dataset_filepath)
    raw_ds.save()
    
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

def get_supplier_infomation(metadata_dict):
    # Connect supplier with material if there is a supplier available
    if "supplier" in metadata_dict:
        
        available_suppliers_openbis = [supplier.props.all()['$name'] for supplier in session.get_samples(type = "SUPPLIER")]
        
        # Create Supplier object or get it from openBIS (in case it is already there)
        if metadata_dict["supplier"]["$name"] not in available_suppliers_openbis:
            supplier_object = create_object_openbis(
                session, 
                "SUPPLIER", 
                "/INSTITUTIONS/SUPPLIERS/SUPPLIER_COLLECTION", 
                metadata_dict["supplier"],
                parents = []
            )
        else:
            supplier_object = get_openbis_object_by_property(session,
                                                            "SUPPLIER", 
                                                            metadata_dict["supplier"]["$name"]
            )
        
        # Remove supplier information from the material metadata because it should be a connection to a SUPPLIER object and not a parameter
        metadata_dict.pop("supplier")
        
        # Append supplier object to the property list
        metadata_dict['has_supplier'] = supplier_object.permId
        
    return metadata_dict

def get_chemist_information(metadata_dict):
    # Connect chemist with material if there is a chemist available
    if "chemist" in metadata_dict:
        
        available_chemists_openbis = [chemist.props.all()['$name'] for chemist in session.get_samples(type = "CHEMIST")]
        
        # Create Chemist object or get it from openBIS (in case it is already there)
        if metadata_dict["chemist"]["$name"] not in available_chemists_openbis:
            chemist_object = create_object_openbis(
                session, 
                "CHEMIST", 
                "/PERSONS/PERSONS/CHEMIST_COLLECTION", 
                metadata_dict["chemist"],
                parents = []
            )
        else:
            chemist_object = get_openbis_object_by_property(session, 
                                                            "CHEMIST", 
                                                            metadata_dict["chemist"]["$name"]
            )
        
        # Remove chemist information from the material metadata because it should be a connection to a CHEMIST object and not a parameter
        metadata_dict.pop("chemist")
        
        # Append supplier object to the property list
        metadata_dict['has_chemist'] = chemist_object.permId
    
    return metadata_dict

def get_molecule_properties(metadata_dict, molecule_identifier, item):
    if is_nan(item["NAMES, FORMULA"]) == False:
        if item["NAMES, FORMULA"].replace(" ", "") not in ["?", ""]: # There are names with
            metadata_dict["$name"] = item["NAMES, FORMULA"]
    else:
        metadata_dict["$name"] = molecule_identifier

    if is_nan(item["DETAILS, COMMENTS"]) == False:
        metadata_dict["comments"] = item["DETAILS, COMMENTS"]

    if is_nan(item["DATE"]) == False:
            processed_date = process_receive_date(item["DATE"])
            if processed_date is not None:
                metadata_dict["receive_date"] = processed_date

    if is_nan(item["STORAGE"]) == False:
        metadata_dict["other_storage_condition"] = True
        metadata_dict["other_storage_condition_specification"] = item["STORAGE"]

    if is_nan(item["ORIGIN"]) == False:
        if item["ORIGIN"]!="?" and item["ORIGIN"].strip()!="":
            metadata_dict["supplier"] = {"$name": item["ORIGIN"]}

    if is_nan(item["SYNTHESIZED BY"]) == False:
        if item["SYNTHESIZED BY"]!="?" and item["SYNTHESIZED BY"].strip()!="":
            metadata_dict["chemist"] = {"$name": item["SYNTHESIZED BY"]}
    
    return metadata_dict
    
def get_substances(filename: str) -> pd.DataFrame: 
    df = pd.read_excel(filename, header = 7)
    return df

def create_molecule_image(chem_mol):
    img = Draw.MolToImage(chem_mol)
    tmpdir = tempfile.mkdtemp()
    struc_filepath = tmpdir + "/" + 'struc.png'
    img.save(struc_filepath)
    
    return struc_filepath

def get_molecule_identifiers(molecule_filepath):
    molecule_identifier = molecule_filepath.split("/")[-1].split(".")[0]
    molecule_number = int(molecule_identifier[:-1])
    return molecule_identifier, molecule_number

def process_receive_date(receive_date):
    if isinstance(receive_date, str):
        processed_receive_date = receive_date.replace(".","/")
        processed_receive_date_split = processed_receive_date.split("/")
        if len(processed_receive_date_split) == 3:
            day, month, year = processed_receive_date.split("/")
            day, month, year = int(day), int(month), int(year)
            # The excel document contains different ways of writing the date...
            if month <= 12:
                processed_receive_date = f"{month}/{day}/{year}"
            else:
                processed_receive_date = f"{day}/{month}/{year}"
        else:
            processed_receive_date = None
    else:
        processed_receive_date = receive_date.strftime('%m/%d/%Y')
    
    return processed_receive_date

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
    parser.add_argument('-cdxml', '--molecules_structures_path', type=str, help='Folder with CDXML files of the molecules', default = None)

    # Read the arguments
    args = parser.parse_args()
    openbis_url = args.openbis_url
    openbis_user = args.openbis_user
    openbis_pw = args.openbis_pw
    material_type = args.material_type
    materials_filepath = args.materials_filepath
    molecules_structures_path = args.molecules_structures_path
    
    if material_type:
        
        # Convert material type to lowercase for an easier processing
        material_type = material_type.lower()
        
        if material_type in MATERIALS_TYPES_LIST:
            
            # Verify if the inventory filepath was entered
            if materials_filepath:
                
                # Flag to check if all the information needed for entering molecules was inserted
                all_input_correct = True
                
                # Connect to openBIS
                session = log_in(bisurl=openbis_url, bisuser=openbis_user, bispasswd=openbis_pw)

                # 2D-layer materials
                if material_type == "2d-layer material":
                    df = pd.read_excel(materials_filepath, header = 2)
                    material_openbis_object_type = "2D_LAYERED_MATERIAL"
                    available_materials_openbis = [material.props.all()["$name"] for material in session.get_samples(type = material_openbis_object_type)]
                
                # Crystals
                elif material_type == "crystal":
                    df = pd.read_excel(materials_filepath, sheet_name = "Drawer")
                    material_openbis_object_type = "CRYSTAL"
                    available_materials_openbis = [material.props.all()["$name"] for material in session.get_samples(type = material_openbis_object_type)]
                
                # Molecules
                elif material_type == "molecule":
                    if molecules_structures_path:
                        df = get_substances(materials_filepath)
                        molecules_filepaths = get_full_filepath(molecules_structures_path)
                        material_openbis_object_type = "MOLECULE"
                        available_materials_openbis = [f"{molecule.props.all()['empa_number']}{molecule.props.all()['batch']}" for molecule in session.get_samples(type = material_openbis_object_type)]
                    else:
                        all_input_correct = False
                
                # In case all inputs were correctly entered, lets upload all the data to openBIS
                if all_input_correct:   
                    
                    for column, item in df.iterrows():
                        if material_type == "2d-layer material":
                            material_name = item["Acronym"]
                        elif material_type == "crystal":
                            material_name = item["Elog Name"]
                        elif material_type == "molecule":
                            molecule_cdxml_found = False
                            
                            for molecule_filepath in molecules_filepaths:
                                # This gets the identifier from the cdxml file (it always returns batch a)
                                molecule_identifier, molecule_number = get_molecule_identifiers(molecule_filepath) # Molecule EMPA ID and EMPA number
                                
                                if is_nan(item["SERIES NUMBER"]) == False:
                                    if int(item["SERIES NUMBER"][0:3]) == molecule_number:
                                        # Correct the molecule batch letter
                                        molecule_batch = item["SERIES NUMBER"][3]
                                        molecule_identifier = f"{molecule_number}{molecule_batch}"
                                        molecule_cdxml_found = True
                                        break  
                            
                            # Identifier is used as a name to generalise (Crystal and 2D-layer material use name to verify if they are already in openBIS.)
                            material_name = molecule_identifier
                        
                        if material_name not in available_materials_openbis:
                            # For the moment I am going to skip the materials with no name. However, this should be solved in the future because every object must have a name.
                            if is_nan(material_name) == False:
                                # Molecule object is much more complex than the other materials. Therefore, it is separated.
                                if material_type == "molecule":
                                    molecule_metadata_dict = {}
                                    # It does not make sense to continue processing the molecule if the cdxml file does not exist
                                    if molecule_cdxml_found:
                                        cdxml_molecule = read_file(molecule_filepath)
                                        molecules = rdkit.Chem.MolsFromCDXML(cdxml_molecule)
                                        
                                        if len(molecules) == 1:
                                            mol = molecules[0] # Get first molecule
                                            mol_chemical_formula = rdMolDescriptors.CalcMolFormula(mol) # Sum Formula
                                            mol_smiles = rdkit.Chem.MolToSmiles(mol) # Canonical Smiles
                                            
                                            material_metadata_dict = {
                                                "empa_number": molecule_number,
                                                "batch": molecule_batch,
                                                "sum_formula": mol_chemical_formula,
                                                "smiles": mol_smiles
                                            }
                                            
                                            material_metadata_dict = get_molecule_properties(material_metadata_dict, molecule_identifier, item)
                                            material_metadata_dict = get_supplier_infomation(material_metadata_dict)
                                            material_metadata_dict = get_chemist_information(material_metadata_dict)
                                            
                                            # Create molecule object
                                            molecule_object = create_object_openbis(
                                                session, 
                                                "MOLECULE", 
                                                f"/MATERIALS/{material_openbis_object_type}S/{material_openbis_object_type}_COLLECTION", 
                                                material_metadata_dict,
                                                parents = []
                                            )
                                            
                                            # Generate molecule image
                                            chem_mol = rdkit.Chem.MolFromSmiles(mol_smiles)
                                            
                                            if chem_mol is not None:
                                                AllChem.Compute2DCoords(chem_mol) # Add coords to the atoms in the molecule
                                                struc_filepath = create_molecule_image(chem_mol)
                                                
                                                # Send molecule image to openBIS
                                                create_dataset_openbis(
                                                    session, 
                                                    'ELN_PREVIEW',
                                                    molecule_object, 
                                                    struc_filepath
                                                )
                                                
                                            else:
                                                print(f"Cannot generate molecule image for {molecule_identifier}")
                                                
                                            # Send molecule cdxml to openBIS
                                            create_dataset_openbis(session, 'RAW_DATA', molecule_object, molecule_filepath)
                                        
                                        elif len(molecules) > 1:
                                            print(f"There are more than one molecule in the file: {molecule_filepath}")  
                                        else:
                                            print(f"There are no molecules in the file: {molecule_filepath}")
                                    
                                    else:
                                        print(f"The CDXML of the molecule {item['SERIES NUMBER']} was not found!")
                                            
                                else:
                                    material_metadata_dict = {"$name": material_name}
                                    
                                    if material_type == "2d-layer material":
                                        material_metadata_dict = get_2d_layer_properties(material_metadata_dict, item)
                                    elif material_type == "crystal":
                                        material_metadata_dict = get_crystal_properties(material_metadata_dict, item)
                                        
                                    material_metadata_dict = get_supplier_infomation(material_metadata_dict)
                                    
                                    material_object = create_object_openbis(
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
                    print("Folder with molecule structures was not introduced.")
                    print(f'Usage: python3 import_materials_to_openbis.py -o <OPENBIS_URL> -u <OPENBIS_USER> -pw <OPENBIS_PW> -type <MATERIAL_TYPE> -xls <MATERIAL_METADATA_EXCEL_FILEPATH> -cdxml <MOLECULE_CDXML_FILES_FOLDERPATH>')
        else:
            print("Unrecognised material types. Allowed material types: Crystal, 2D-Layer Material, ...")
            print(f'Usage: python3 import_materials_to_openbis.py -o <OPENBIS_URL> -u <OPENBIS_USER> -pw <OPENBIS_PW> -type <MATERIAL_TYPE> -xls <MATERIAL_METADATA_EXCEL_FILEPATH> -cdxml <MOLECULE_CDXML_FILES_FOLDERPATH>')
    else:
        print("Specify the material type. Ex.: Crystal, Molecule, 2D-Layer Material, ...")
        print(f'Usage: python3 import_materials_to_openbis.py -o <OPENBIS_URL> -u <OPENBIS_USER> -pw <OPENBIS_PW> -type <MATERIAL_TYPE> -xls <MATERIAL_METADATA_EXCEL_FILEPATH> -cdxml <MOLECULE_CDXML_FILES_FOLDERPATH>')