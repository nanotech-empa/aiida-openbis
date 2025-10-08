from pybis import Openbis
import argparse
import os
import pandas as pd
import json


def log_in(bisurl="openbis", bisuser="admin", bispasswd="changeit"):
    """Function to login to openBIS."""
    if Openbis(bisurl, verify_certificates=False).is_token_valid():
        session = Openbis(bisurl, verify_certificates=False)
    else:
        Openbis(bisurl, verify_certificates=False).login(
            bisuser, bispasswd, save_token=True
        )
        session = Openbis(bisurl, verify_certificates=False)
    return session


def create_json(json_content, filename):
    with open(filename, "w") as file:
        json.dump(json_content, file, indent=4)


def upload_instruments_to_openbis(session, folderpath, save_config_folderpath):
    data_filepath = os.path.join(folderpath, "data.xlsx")
    instruments_data = pd.read_excel(
        data_filepath, header=0, dtype=str, sheet_name="Instruments"
    )
    instruments_data = instruments_data.fillna("")
    components_data = pd.read_excel(
        data_filepath, header=0, dtype=str, sheet_name="Components"
    )
    components_data = components_data.fillna("")

    instruments_metadata = {}
    for instrument_idx, row in instruments_data.iterrows():
        components_indexes = row["components"].split(",")
        components_identifiers = []
        components_metadata = {}
        for component_idx in components_indexes:
            component_data = components_data[components_data["index"] == component_idx]
            component_name = component_data["name"].values[0]
            component_actions = component_data["actions"].values[0]
            component_actions = (
                component_actions.replace(" ", "").split(",")
                if component_actions
                else []
            )
            component_observables = component_data["observables"].values[0]
            component_observables = (
                component_observables.replace(" ", "").split(",")
                if component_observables
                else []
            )
            component_settings = component_data["settings"].values[0]
            component_settings = (
                component_settings.replace(" ", "").split(",")
                if component_settings
                else []
            )

            component_props = {"name": component_name}

            component_object = session.new_sample(
                type="COMPONENT",
                experiment="/EQUIPMENT/ILOG/COMPONENT_COLLECTION",
                props=component_props,
            )
            component_object.save()
            component_identifier = component_object.permId
            components_identifiers.append(component_identifier)

            components_metadata[component_identifier] = {
                "name": component_name,
                "actions": component_actions,
                "observables": component_observables,
                "settings": component_settings,
            }

        instrument_name = row["name"]
        instrument_props = {
            "name": instrument_name,
            "components": components_identifiers,
        }

        instrument_object = session.new_sample(
            type="INSTRUMENT",
            experiment="/EQUIPMENT/ILOG/INSTRUMENT_COLLECTION",
            props=instrument_props,
        )
        instrument_object.save()

        instrument_identifier = instrument_object.permId

        instruments_metadata[instrument_identifier] = {
            "name": instrument_name,
            "components": components_metadata,
        }

    create_json(
        instruments_metadata, f"{save_config_folderpath}/instruments_config.json"
    )


# Set up openBIS environment
if __name__ == "__main__":
    # Read command line arguments

    parser = argparse.ArgumentParser(
        description="Upload instruments into openBIS and create instrument config file for webapp interface."
    )

    # Define the arguments with flags
    parser.add_argument(
        "-o",
        "--openbis_url",
        type=str,
        help="OpenBIS URL",
        default="https://local.openbis.ch/openbis",
    )
    parser.add_argument(
        "-u", "--openbis_user", type=str, help="OpenBIS User", default="admin"
    )
    parser.add_argument(
        "-pw", "--openbis_pw", type=str, help="OpenBIS Password", default="123456789"
    )
    parser.add_argument(
        "-dir",
        "--instruments_folderpath",
        type=str,
        help="Directory with instrument information",
        default=None,
    )
    parser.add_argument(
        "-save_dir",
        "--save_config_folderpath",
        type=str,
        help="Directory to save instruments config file",
        default="/home/jovyan/aiida-openbis/Notebooks/connection_to_openbis",
    )

    # Read the arguments
    args = parser.parse_args()
    openbis_url = args.openbis_url
    openbis_user = args.openbis_user
    openbis_pw = args.openbis_pw
    instruments_folderpath = args.instruments_folderpath
    save_instruments_config_folderpath = args.save_config_folderpath

    # Open openBIS session
    openbis_session = log_in(openbis_url, openbis_user, openbis_pw)

    # Upload instruments to openBIS and make instruments_config file for the webapp
    upload_instruments_to_openbis(
        openbis_session, instruments_folderpath, save_instruments_config_folderpath
    )
