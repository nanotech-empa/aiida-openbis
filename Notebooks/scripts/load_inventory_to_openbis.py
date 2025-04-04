from pybis import Openbis
import argparse
import os
import pandas as pd
import rdkit
import json
from rdkit.Chem import AllChem, Draw, rdMolDescriptors
import tempfile

def read_file(filepath: str) -> str:
    with open(filepath, "rb") as f:
        return f.read()

def create_json(json_content, filename):
    with open(filename, 'w') as file:
        json.dump(json_content, file, indent=4)

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

def upload_molecules_to_openbis(session, datasets_folderpath):
    data_filepath = os.path.join(datasets_folderpath, "inventory_data.xlsx")
    structures_folderpath = os.path.join(datasets_folderpath, "molecules_dataset", "structures")
    structures_workshop_folderpath = os.path.join(datasets_folderpath, "molecules_dataset", "structures_workshop")
    structures_filepaths = get_full_filepath(structures_folderpath) + get_full_filepath(structures_workshop_folderpath)
    data = pd.read_excel(data_filepath, header = 0, dtype=str, sheet_name = "Molecule")
    data = data.fillna('')
    
    for _, row in data.iterrows():
        molecule_props = {}
        molecule_concept_props = {}
        
        # ------ CREATE MOLECULE CONCEPT OBJECT ------
        molecule_series_number = row["empa_number"]
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
        molecule_batch = row["batch"]
        molecule_vial = row["vial"]
        if row["name"]:
            molecule_name = row["name"]
        else:
            molecule_name = f"{molecule_series_number}{molecule_batch}{molecule_vial}"
            
        molecule_receive_date = convert_date(row["receive_date"], "%m/%d/%Y")
        
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

def upload_instruments_to_openbis(session, folderpath):
    data_filepath = os.path.join(folderpath, "inventory_data.xlsx")
    instruments_data = pd.read_excel(data_filepath, header = 0, dtype=str, sheet_name = "Instrument")
    instruments_data = instruments_data.fillna('')
    components_data = pd.read_excel(data_filepath, header = 0, dtype=str, sheet_name = "Component")
    components_data = components_data.fillna('')
    
    instruments_metadata = {}
    for instrument_idx, row in instruments_data.iterrows():
        if row["components"]:
            components_indexes = row["components"].split(",")
        else:
            continue
        
        components_identifiers = []
        components_metadata = {}
        for component_idx in components_indexes:
            component_data = components_data[components_data["index"] == component_idx]
            component_name = component_data["name"].values[0]
            component_actions = component_data["actions"].values[0]
            component_actions = component_actions.replace(" ", "").split(",") if component_actions else []
            component_observables = component_data["observables"].values[0]
            component_observables = component_observables.replace(" ", "").split(",") if component_observables else []
            component_settings = component_data["settings"].values[0]
            component_settings = component_settings.replace(" ", "").split(",") if component_settings else []
            
            component_props = {
                "name": component_name
            }
            
            component_object = session.new_sample(
                type = "COMPONENT", 
                experiment = "/EQUIPMENT/ILOG/COMPONENT_COLLECTION", 
                props = component_props
            )
            component_object.save()
            component_identifier = component_object.permId
            components_identifiers.append(component_identifier)
            
            components_metadata[component_identifier] = {
                "name": component_name,
                "actions": component_actions,
                "observables": component_observables,
                "settings": component_settings
            }
        
        instrument_name = row["name"]
        instrument_props = {
            "name": instrument_name,
            "components": components_identifiers
        }
        
        instrument_object = session.new_sample(
            type = "INSTRUMENT", 
            experiment = "/EQUIPMENT/ILOG/INSTRUMENT_COLLECTION", 
            props = instrument_props
        )
        instrument_object.save()
        
        instrument_identifier = instrument_object.permId
        
        instruments_metadata[instrument_identifier] = {
            "name": instrument_name,
            "components": components_metadata
        }
        
    create_json(instruments_metadata, f"/home/jovyan/aiida-openbis/Notebooks/connection_to_openbis/instruments_config.json")

def upload_inventory_to_openbis(session, folderpath):
    data_filepath = os.path.join(folderpath, "inventory_data.xlsx")
    grants_data = pd.read_excel(data_filepath, header = 0, dtype=str, sheet_name = "Grant")
    grants_data = grants_data.fillna('')
    
    institutions_data = pd.read_excel(data_filepath, header = 0, dtype=str, sheet_name = "Institution")
    institutions_data = institutions_data.fillna('')
    
    persons_data = pd.read_excel(data_filepath, header = 0, dtype=str, sheet_name = "Person")
    persons_data = persons_data.fillna('')
    
    protocols_data = pd.read_excel(data_filepath, header = 0, dtype=str, sheet_name = "Protocol")
    protocols_data = protocols_data.fillna('')
    
    chemical_concepts_data = pd.read_excel(data_filepath, header = 0, dtype=str, sheet_name = "ChemicalConcept")
    chemical_concepts_data = chemical_concepts_data.fillna('')
    
    chemicals_data = pd.read_excel(data_filepath, header = 0, dtype=str, sheet_name = "Chemical")
    chemicals_data = chemicals_data.fillna('')
    
    crystal_concepts_data = pd.read_excel(data_filepath, header = 0, dtype=str, sheet_name = "CrystalConcept")
    crystal_concepts_data = crystal_concepts_data.fillna('')
    
    crystals_data = pd.read_excel(data_filepath, header = 0, dtype=str, sheet_name = "Crystal")
    crystals_data = crystals_data.fillna('')
    
    reaction_product_concepts_data = pd.read_excel(data_filepath, header = 0, dtype=str, sheet_name = "ReactionProductConcept")
    reaction_product_concepts_data = reaction_product_concepts_data.fillna('')
    
    reaction_products_data = pd.read_excel(data_filepath, header = 0, dtype=str, sheet_name = "ReactionProduct")
    reaction_products_data = reaction_products_data.fillna('')
    
    wafer_substrate_concepts_data = pd.read_excel(data_filepath, header = 0, dtype=str, sheet_name = "WaferSubstrateConcept")
    wafer_substrate_concepts_data = wafer_substrate_concepts_data.fillna('')
    
    wafer_substrates_data = pd.read_excel(data_filepath, header = 0, dtype=str, sheet_name = "WaferSubstrate")
    wafer_substrates_data = wafer_substrates_data.fillna('')
    
    twod_material_concepts_data = pd.read_excel(data_filepath, header = 0, dtype=str, sheet_name = "TwoDLayerMaterialConcept")
    twod_material_concepts_data = twod_material_concepts_data.fillna('')
    
    twod_materials_data = pd.read_excel(data_filepath, header = 0, dtype=str, sheet_name = "TwoDLayerMaterial")
    twod_materials_data = twod_materials_data.fillna('')
    
    # Grants
    for idx, row in grants_data.iterrows():
        props = {
            "name": row["name"],
            "funder": row["funder"],
        }
        object = session.new_sample(
            type = "GRANT", 
            experiment = "/ADMINISTRATIVE/FUNDING/GRANT_COLLECTION", 
            props = props
        )
        object.save()
    
    # Institutions
    institutions_objects = []
    for idx, row in institutions_data.iterrows():
        props = {
            "name": row["name"],
        }
        object = session.new_sample(
            type = "INSTITUTION", 
            experiment = "/ADMINISTRATIVE/INSTITUTIONS/ACADEMIC_INSTITUTION_COLLECTION", 
            props = props
        )
        object.save()
        
        institutions_objects.append(object)
        
    # Persons
    for idx, row in persons_data.iterrows():
        institution_index = int(row["institution"])
        props = {
            "name": row["name"],
            "institutions": [institutions_objects[institution_index - 1].permId] # 0-index correction
        }
        object = session.new_sample(
            type = "PERSON", 
            experiment = "/ADMINISTRATIVE/PERSONS/LAB_205_PERSON_COLLECTION", 
            props = props
        )
        object.save()
    
    # Protocols
    for idx, row in protocols_data.iterrows():
        props = {
            "name": row["name"],
            "short_name": row["short_name"],
            "protocol_template": row["protocol_template"]
        }
        object = session.new_sample(
            type = "PROTOCOL", 
            experiment = "/METHODS/PROTOCOLS/PROTOCOL_COLLECTION", 
            props = props
        )
        object.save()
    
    # Chemical Concepts
    chemical_concept_objects = []
    for idx, row in chemical_concepts_data.iterrows():
        props = {
            "name": row["name"]
        }
        object = session.new_sample(
            type = "CHEMICAL_CONCEPT", 
            experiment = "/MATERIALS/RAW_MATERIALS/CHEMICAL_COLLECTION", 
            props = props
        )
        object.save()
        chemical_concept_objects.append(object)
    
    # Chemicals
    for idx, row in chemicals_data.iterrows():
        chemical_concept = int(row["chemical_concept"])
        props = {
            "name": row["name"],
            "chemical_concept": chemical_concept_objects[chemical_concept - 1].permId
        }
        object = session.new_sample(
            type = "CHEMICAL", 
            experiment = "/MATERIALS/RAW_MATERIALS/CHEMICAL_COLLECTION", 
            props = props
        )
        object.save()
    
    # Crystal Concepts
    crystal_concept_objects = []
    for idx, row in crystal_concepts_data.iterrows():
        props = {
            "name": row["name"],
            "material": row["material"],
            "face": row["face"]
        }
        object = session.new_sample(
            type = "CRYSTAL_CONCEPT", 
            experiment = "/MATERIALS/CRYSTALS/AU_COLLECTION", # TODO: Put in the right collection
            props = props
        )
        object.save()
        crystal_concept_objects.append(object)
    
    # Crystals
    for idx, row in crystals_data.iterrows():
        crystal_concept = int(row["crystal_concept"])
        props = {
            "name": row["name"],
            "crystal_concept": crystal_concept_objects[crystal_concept - 1].permId
        }
        object = session.new_sample(
            type = "CRYSTAL", 
            experiment = "/MATERIALS/CRYSTALS/AU_COLLECTION", 
            props = props
        )
        object.save()
        
# Set up openBIS environment
if __name__ == "__main__":
    # Read command line arguments
    
    parser = argparse.ArgumentParser(description = 'Setup the openBIS database (create objects types, collections, etc.).')

    # Define the arguments with flags
    parser.add_argument('-o', '--openbis_url', type=str, help='OpenBIS URL', default = 'https://local.openbis.ch/openbis')
    parser.add_argument('-u', '--openbis_user', type=str, help='OpenBIS User', default = 'admin')
    parser.add_argument('-pw', '--openbis_pw', type=str, help='OpenBIS Password', default = '123456789')
    parser.add_argument('-dir', '--datasets_folderpath', type=str, help='Directory with inventory information', default = None)

    # Read the arguments
    args = parser.parse_args()
    openbis_url = args.openbis_url
    openbis_user = args.openbis_user
    openbis_pw = args.openbis_pw
    datasets_folderpath = args.datasets_folderpath
    
    # Open openBIS session
    openbis_session = log_in(openbis_url, openbis_user, openbis_pw)
    
    # Upload molecules to openBIS
    upload_molecules_to_openbis(openbis_session, datasets_folderpath)
    
    # Upload instruments to openBIS
    upload_instruments_to_openbis(openbis_session, datasets_folderpath)
    
    # Upload inventory to openBIS
    upload_inventory_to_openbis(openbis_session, datasets_folderpath)