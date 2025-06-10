from pybis import Openbis
import argparse
import os
import pandas as pd
import rdkit
import json
from rdkit.Chem import AllChem, Draw, rdMolDescriptors
import tempfile

def log_in(bisurl='openbis', bisuser='admin', bispasswd='changeit'):
    """Function to login to openBIS."""
    if Openbis(bisurl, verify_certificates=False).is_token_valid():
        session = Openbis(bisurl, verify_certificates=False)
    else:
        Openbis(bisurl, verify_certificates=False).login(bisuser, bispasswd, save_token=True)
        session = Openbis(bisurl, verify_certificates=False)
    return session

def upload_inventories(session, dataset_folderpath):
    datasets_map = {
        "people.xlsx": ["PERSON", "/ADMINISTRATIVE/PEOPLE/LAB_205_PERSON_COLLECTION"],
        "institutions.xlsx": ["INSTITUTION", "/ADMINISTRATIVE/INSTITUTIONS/ACADEMIC_INSTITUTION_COLLECTION"]
    }
    
    full_filepaths = [os.path.join(dataset_folderpath, file) for file in os.listdir(dataset_folderpath)]

# Set up openBIS environment
if __name__ == "__main__":
    # Read command line arguments
    
    parser = argparse.ArgumentParser(description = 'Setup the openBIS database (create objects types, collections, etc.).')

    # Define the arguments with flags
    parser.add_argument('-o', '--openbis_url', type=str, help='OpenBIS URL', default = 'https://local.openbis.ch/openbis')
    parser.add_argument('-u', '--openbis_user', type=str, help='OpenBIS User', default = 'admin')
    parser.add_argument('-pw', '--openbis_pw', type=str, help='OpenBIS Password', default = '123456789')
    parser.add_argument('-token', '--openbis_token', type=str, help='OpenBIS Token', default = None)
    parser.add_argument('-dir', '--datasets_folderpath', type=str, help='Directory with inventory information', default = None)

    # Read the arguments
    args = parser.parse_args()
    openbis_url = args.openbis_url
    openbis_user = args.openbis_user
    openbis_pw = args.openbis_pw
    openbis_token = args.openbis_token
    datasets_folderpath = args.datasets_folderpath
    
    # Open openBIS session
    openbis_session = log_in(openbis_url, openbis_user, openbis_pw)
    
    # Upload inventories to openBIS
    upload_inventories(openbis_session, datasets_folderpath)