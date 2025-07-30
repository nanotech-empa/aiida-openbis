from pybis import Openbis
import argparse
import os
import pandas as pd
import rdkit
from rdkit.Chem import AllChem, Draw, rdMolDescriptors
import tempfile

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
    return [f"{folderpath}/{file}" for file in os.listdir(folderpath)]

def find_cdxml_file(filepaths, number):
    # Convert the number to string and check both padded and non-padded formats
    num_str = str(number)
    num_padded = num_str.zfill(3)  # Pads with zeros to make it 3 digits
    
    # Search for the file with exact filename match (before .cdxml)
    for path in filepaths:
        filename = os.path.splitext(os.path.basename(path))[0]  # Get the filename without extension
        if filename == num_str or filename == num_padded:
            return path
    
    # If no match is found, return None
    return None

def get_openbis_object_by_name(session, object_type, property_value):
    for object in session.get_samples(type = object_type):
        if object.props.all()['name'] == property_value:
            return object
        
def convert_date(date_str, format):
    try:
        return pd.to_datetime(date_str).strftime("%m/%d/%Y")
    except ValueError:
        return None
    
def create_molecule_image(chem_mol):
    img = Draw.MolToImage(chem_mol)
    tmpdir = tempfile.mkdtemp()
    struc_filepath = tmpdir + "/" + 'struc.png'
    img.save(struc_filepath)
    return struc_filepath

def upload_molecules_to_openbis(session, molecules_folderpath):
    data_filepath = os.path.join(molecules_folderpath, "data.xlsx")
    structures_folderpath = os.path.join(molecules_folderpath, "structures")
    structures_workshop_folderpath = os.path.join(molecules_folderpath, "structures_workshop")
    structures_filepaths = get_full_filepath(structures_folderpath) + get_full_filepath(structures_workshop_folderpath)
    data = pd.read_excel(data_filepath, header = 0, dtype=str)
    data = data.fillna('')
    
    for _, row in data.iterrows():
        molecule_props = {}
        molecule_concept_props = {}
        
        # ------ CREATE MOLECULE CONCEPT OBJECT ------
        molecule_series_number = row["Series Number"]
        molecule_concept_props["name"] = molecule_series_number
        
        # Find molecule concept object in openBIS
        molecule_concept_object = get_openbis_object_by_name(session, "MOLECULE_CONCEPT", molecule_series_number)
        
        if molecule_concept_object:
            print("Molecule concept is already in openBIS.")
        else:
            # Get molecule concept properties
            molecule_cdxml_filepath = find_cdxml_file(structures_filepaths, molecule_series_number)
            molecule_image_filepath = ""
            
            if molecule_cdxml_filepath:
                cdxml_molecule = read_file(molecule_cdxml_filepath)
                molecules = rdkit.Chem.MolsFromCDXML(cdxml_molecule)
            
                if len(molecules) == 1:
                    mol = molecules[0] # Get first molecule
                    mol_chemical_formula = rdMolDescriptors.CalcMolFormula(mol) # Sum Formula
                    mol_smiles = rdkit.Chem.MolToSmiles(mol) # Canonical Smiles
                    molecule_concept_props["smiles"] = mol_smiles
                    molecule_concept_props["sum_formula"] = mol_chemical_formula
                    
                    # Generate molecule image
                    molecule_image_filepath = create_molecule_image(mol)
                    
                elif len(molecules) > 1:
                    print(f"There are more than one molecule in the file: {molecule_cdxml_filepath}")  
                else:
                    print(f"There are no molecules in the file: {molecule_cdxml_filepath}")
            
            molecule_concept_object = session.new_sample(
                type = "MOLECULE_CONCEPT", 
                experiment = "/MATERIALS/MOLECULES/PRECURSOR_COLLECTION", 
                props = molecule_concept_props
            )
            molecule_concept_object.save()
            
            if molecule_cdxml_filepath:
                session.new_dataset(type = "RAW_DATA", object = molecule_concept_object, file = molecule_cdxml_filepath).save()
            if molecule_image_filepath:
                session.new_dataset(type = "ELN_PREVIEW", object = molecule_concept_object, file = molecule_image_filepath).save()
        
        # ------ CREATE MOLECULE OBJECT ------
        molecule_batch = row["Batch"]
        molecule_vial = row["Vial"]
        if row["Name"]:
            molecule_name = row["Name"]
        else:
            molecule_name = f"{molecule_series_number}{molecule_batch}{molecule_vial}"
            
        molecule_receive_date = convert_date(row["Date"], "%m/%d/%Y")
        
        # Find molecule object in openBIS
        molecule_object = get_openbis_object_by_name(session, "MOLECULE", molecule_name)
        
        if molecule_object:
            print("Molecule is already in openBIS.")
            continue
        else:
            # Molecule properties
            molecule_props["name"] = molecule_name
            molecule_props["empa_number"] = int(molecule_series_number)
            molecule_props["molecule_concept"] = molecule_concept_object.permId
            if molecule_batch:
                molecule_props["batch"] = molecule_batch
            if molecule_vial:
                molecule_props["vial"] = molecule_vial
            if molecule_receive_date:
                molecule_props["receive_date"] = molecule_receive_date
            
            molecule_object = session.new_sample(
                type = "MOLECULE", 
                experiment = "/MATERIALS/MOLECULES/PRECURSOR_COLLECTION", 
                props = molecule_props
            )
            molecule_object.save()

# Set up openBIS environment
if __name__ == "__main__":
    # Read command line arguments
    
    parser = argparse.ArgumentParser(description = 'Setup the openBIS database (create objects types, collections, etc.).')

    # Define the arguments with flags
    parser.add_argument('-o', '--openbis_url', type=str, help='OpenBIS URL', default = 'https://local.openbis.ch/openbis')
    parser.add_argument('-u', '--openbis_user', type=str, help='OpenBIS User', default = 'admin')
    parser.add_argument('-pw', '--openbis_pw', type=str, help='OpenBIS Password', default = '123456789')
    parser.add_argument('-dir', '--molecules_folderpath', type=str, help='Directory with molecule information', default = None)

    # Read the arguments
    args = parser.parse_args()
    openbis_url = args.openbis_url
    openbis_user = args.openbis_user
    openbis_pw = args.openbis_pw
    molecules_folderpath = args.molecules_folderpath
    
    # Open openBIS session
    openbis_session = log_in(openbis_url, openbis_user, openbis_pw)
    
    # Upload molecules to openBIS
    upload_molecules_to_openbis(openbis_session, molecules_folderpath)