# Load libraries
from aiida_openbis.utils import bisutils
import utils
import os

# Connect to openBIS
session = bisutils.log_in(bisurl="openbis", bisuser="admin", bispasswd="changeit")

# Load property types and object types from YAML file

objects_folder = '/home/jovyan/work/aiida-openbis/Python_Scripts/Objects/'
config_filenames = os.listdir(objects_folder)

for config_filename in config_filenames:

    if '.yml' in config_filename:
        mode = 'r'
        filepath = f'/home/jovyan/work/aiida-openbis/Python_Scripts/Objects/{config_filename}'
        items = utils.load_yml_file(filepath, mode)

        # Create new vocabulary
        if 'vocabularies' in items:
            utils.create_new_vocabulary(session, items)

        # Create new property types
        if 'property_types' in items:
            utils.create_new_property_types(session, items)

        # Create new object types
        if 'object_types' in items:
            utils.create_new_object_types(session, items)