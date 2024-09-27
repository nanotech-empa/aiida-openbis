# Import libraries
import os
import re
import pandas as pd
import numpy as np
import rdkit
import tempfile
import datetime
from rdkit.Chem import rdMolDescriptors
from pybis import Openbis
from rdkit.Chem import AllChem, Draw
import argparse
import warnings
warnings.filterwarnings("ignore")

# Functions

def is_nan(value):
    return value != value

def get_substances(filename: str) -> pd.DataFrame: 
    df = pd.read_excel(filename, header = 7)
    return df

def get_full_filepath(folderpath: str) -> list[str]:
    return [f"{folderpath}/{file}" for file in os.listdir(folderpath) if file.endswith(".cdxml")]

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

def create_dataset_openbis(session, dataset_type, object_id, dataset_filepath):
    raw_ds = session.new_dataset(type = dataset_type, object = object_id, file = dataset_filepath)
    raw_ds.save()
    
def create_molecule_image(chem_mol):
    img = Draw.MolToImage(chem_mol)
    tmpdir = tempfile.mkdtemp()
    struc_filepath = tmpdir + "/" + 'struc.png'
    img.save(struc_filepath)
    
    return struc_filepath

def get_openbis_object_by_property(session, object_type, property_value):
    for object in session.get_samples(type = object_type):
        if object.props.all()['$name'] == property_value:
            return object
        
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
    parser.add_argument('-xls', '--substances_path', type=str, help='Excel file with list of substances', default = None)
    parser.add_argument('-cdxml', '--molecules_structures_path', type=str, help='Folder with CDXML files of the molecules', default = None)

    args = parser.parse_args()
    
    openbis_url = args.openbis_url
    openbis_user = args.openbis_user
    openbis_pw = args.openbis_pw
    substances_path = args.substances_path
    molecules_structures_path = args.molecules_structures_path

    if substances_path and molecules_structures_path:
        # Connect to openBIS
        session = log_in(bisurl=openbis_url, bisuser=openbis_user, bispasswd=openbis_pw)

        # Get molecules from cdxml files
        all_substances_info = get_substances(substances_path)
        molecules_filepaths = get_full_filepath(molecules_structures_path)

        molecules_metadata = []

        available_molecules_openbis = [f"{molecule.props.all()['empa_number']}{molecule.props.all()['batch']}" for molecule in session.get_samples(type = "MOLECULE")]

        for substance_idx, molecule_info in all_substances_info.iterrows():
            molecule_cdxml_found = False
            
            for molecule_filepath in molecules_filepaths:
                molecule_identifier, molecule_number = get_molecule_identifiers(molecule_filepath) # Molecule EMPA ID and EMPA number
                
                if is_nan(molecule_info["SERIES NUMBER"]) == False:
                    if int(molecule_info["SERIES NUMBER"][0:3]) == molecule_number:
                        molecule_batch = molecule_info["SERIES NUMBER"][3]
                        molecule_identifier = f"{molecule_number}{molecule_batch}"
                        molecule_cdxml_found = True
                        break
                else:
                    continue # If the series number does not exist, move to the next molecule in the list
            
            if molecule_cdxml_found:
                cdxml_molecule = read_file(molecule_filepath)
                molecules = rdkit.Chem.MolsFromCDXML(cdxml_molecule)
                
                if len(molecules) == 1:
                    mol = molecules[0] # Get first molecule
                    mol_chemical_formula = rdMolDescriptors.CalcMolFormula(mol) # Sum Formula
                    mol_smiles = rdkit.Chem.MolToSmiles(mol) # Canonical Smiles
                    
                    if molecule_identifier not in available_molecules_openbis: # pyBIS returns vocabulary terms in capital letters
                        molecule_metadata_dict = {}
                        
                        if is_nan(molecule_info["NAMES, FORMULA"]) == False:
                            if molecule_info["NAMES, FORMULA"].replace(" ", "") not in ["?", ""]: # There are names with
                                molecule_metadata_dict["$name"] = molecule_info["NAMES, FORMULA"]
                        else:
                            molecule_metadata_dict["$name"] = molecule_identifier

                        if is_nan(molecule_info["DETAILS, COMMENTS"]) == False:
                            molecule_metadata_dict["comments"] = molecule_info["DETAILS, COMMENTS"]

                        if is_nan(molecule_info["DATE"]) == False:
                                processed_date = process_receive_date(molecule_info["DATE"])
                                if processed_date is not None:
                                    molecule_metadata_dict["receive_date"] = processed_date

                        if is_nan(molecule_info["STORAGE"]) == False:
                            molecule_metadata_dict["other_storage_condition"] = True
                            molecule_metadata_dict["other_storage_condition_specification"] = molecule_info["STORAGE"]

                        if is_nan(molecule_info["ORIGIN"]) == False:
                            if molecule_info["ORIGIN"]!="?" and molecule_info["ORIGIN"].strip()!="":
                                molecule_metadata_dict["supplier"] = {"$name": molecule_info["ORIGIN"]}

                        if is_nan(molecule_info["SYNTHESIZED BY"]) == False:
                            if molecule_info["SYNTHESIZED BY"]!="?" and molecule_info["SYNTHESIZED BY"].strip()!="":
                                molecule_metadata_dict["chemist"] = {"$name": molecule_info["SYNTHESIZED BY"]}
                        
                        molecule_metadata_dict["empa_number"] = molecule_number
                        molecule_metadata_dict["batch"] = molecule_batch
                        molecule_metadata_dict["sum_formula"] = mol_chemical_formula
                        molecule_metadata_dict["smiles"] = mol_smiles
                        
                        molecule_parents = []
                        
                        # Connect supplier with molecule if there is a supplier available
                        if "supplier" in molecule_metadata_dict:
                            
                            available_suppliers_openbis = [supplier.props.all()['$name'] for supplier in session.get_samples(type = "SUPPLIER")]
                            
                            # Create Supplier object or get it from openBIS (in case it is already there)
                            if molecule_metadata_dict["supplier"]["$name"] not in available_suppliers_openbis:
                                supplier_object = create_object_openbis(
                                    session, 
                                    "SUPPLIER", 
                                    "/INSTITUTIONS/SUPPLIERS/SUPPLIER_COLLECTION", 
                                    molecule_metadata_dict["supplier"],
                                    parents = []
                                )
                            else:
                                supplier_object = get_openbis_object_by_property(session, 
                                                                                "SUPPLIER", 
                                                                                molecule_metadata_dict["supplier"]["$name"]
                                )
                            
                            # Remove supplier information from the molecule metadata because it should be a connection to a SUPPLIER object and not a parameter
                            molecule_metadata_dict.pop("supplier")
                            
                            # Append supplier object to the property list
                            molecule_metadata_dict['has_supplier'] = supplier_object.permId
                        
                        # Connect chemist with molecule if there is a chemist available
                        if "chemist" in molecule_metadata_dict:
                            
                            available_chemists_openbis = [supplier.props.all()['$name'] for supplier in session.get_samples(type = "CHEMIST")]
                            
                            # Create Chemist object or get it from openBIS (in case it is already there)
                            if molecule_metadata_dict["chemist"]["$name"] not in available_chemists_openbis:
                                chemist_object = create_object_openbis(
                                    session, 
                                    "CHEMIST", 
                                    "/PERSONS/PERSONS/CHEMIST_COLLECTION", 
                                    molecule_metadata_dict["chemist"],
                                    parents = []
                                )
                            else:
                                chemist_object = get_openbis_object_by_property(session, 
                                                                                "CHEMIST", 
                                                                                molecule_metadata_dict["chemist"]["$name"]
                                )
                            
                            # Remove chemist information from the molecule metadata because it should be a connection to a CHEMIST object and not a parameter
                            molecule_metadata_dict.pop("chemist")
                            
                            # Append supplier object to the property list
                            molecule_metadata_dict['has_chemist'] = chemist_object.permId
                        
                        # Create molecule object
                        molecule_object = create_object_openbis(
                            session, 
                            "MOLECULE", 
                            "/MATERIALS/MOLECULES/MOLECULE_COLLECTION", 
                            molecule_metadata_dict,
                            parents = molecule_parents
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
                        create_dataset_openbis(
                            session, 
                            'RAW_DATA', 
                            molecule_object, 
                            molecule_filepath
                        )
                        
                    else:
                        print(f"Molecule {molecule_identifier} already in openBIS.")
                        
                elif len(molecules) > 1:
                    print(f"There are more than one molecule in the file: {molecule_filepath}")  
                else:
                    print(f"There are no molecules in the file: {molecule_filepath}")
            else:
                print(f"The CDXML of the molecule {molecule_info['SERIES NUMBER']} was not found!")
    else:
        print("The path to the excel file with molecule metadata or to the folder with CDXML files were not entered.")
        print(f'Usage: python3 import_molecules_to_openbis.py -o <OPENBIS_URL> -u <OPENBIS_USER> -pw <OPENBIS_PW> -xls <MOLECULE_METADATA_EXCEL_FILEPATH> -cdxml <MOLECULE_CDXML_FILES_FOLDERPATH>')