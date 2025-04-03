from pybis import Openbis
import json
import argparse

def log_in(bisurl='openbis', bisuser='admin', bispasswd='changeit'):
    """Function to login to openBIS."""
    if Openbis(bisurl, verify_certificates=False).is_token_valid():
        session = Openbis(bisurl, verify_certificates=False)
    else:
        Openbis(bisurl, verify_certificates=False).login(bisuser, bispasswd, save_token=True)
        session = Openbis(bisurl, verify_certificates=False)
    return session

def upload_protocols_to_openbis(session):
    protocols_metadata = [
        {
            "name": "Sputter + Anneal",
            "short_name": "SPAN",
            "instrument": "20250401142245758-3054",
            "protocol_template": json.dumps(
                [
                    {
                        "name": "Process 1", 
                        "actions": [
                            {
                                "name": "Sputter 1", 
                                "components": [
                                    "20250401142240887-3025"
                                ], 
                                "action_type": "Sputtering"
                            }
                        ], 
                        "observables": [
                            {
                                "name": "Current", 
                                "component": "20250401142240887-3025", 
                                "observable_type": "Observable"
                            }
                        ]
                    }, 
                    {
                        "name": "Process 2", 
                        "actions": [
                            {
                                "name": "Anneal 1", 
                                "components": [
                                    "20250401142240247-3021"
                                ], 
                                "action_type": "Annealing"
                            }
                        ], 
                        "observables": [
                            {
                                "name": "Current", 
                                "component": "20250401142240247-3021", 
                                "observable_type": "Observable"
                            }
                        ]
                    }
                ]
            )
        }
    ]
    
    for protocol_props in protocols_metadata:
        protocol_object = session.new_sample(
            type = "PROTOCOL", 
            experiment = "/METHODS/PROTOCOLS/PROTOCOL_COLLECTION", 
            props = protocol_props
        )
        protocol_object.save()
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = 'Upload instruments into openBIS and create instrument config file for webapp interface.')

    # Define the arguments with flags
    parser.add_argument('-o', '--openbis_url', type=str, help='OpenBIS URL', default = 'https://local.openbis.ch/openbis')
    parser.add_argument('-u', '--openbis_user', type=str, help='OpenBIS User', default = 'admin')
    parser.add_argument('-pw', '--openbis_pw', type=str, help='OpenBIS Password', default = '123456789')

    # Read the arguments
    args = parser.parse_args()
    openbis_url = args.openbis_url
    openbis_user = args.openbis_user
    openbis_pw = args.openbis_pw
    
    # Open openBIS session
    openbis_session = log_in(openbis_url, openbis_user, openbis_pw)
    
    # Upload instruments to openBIS
    upload_protocols_to_openbis(openbis_session)