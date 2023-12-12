# Load libraries
from aiida_openbis.utils import bisutils
import utils
import os

# Connect to openBIS
session = bisutils.log_in(bisurl="openbis", bisuser="admin", bispasswd="changeit")

# Load property types and object types from YAML file

objects_folder = '/home/jovyan/work/aiida-openbis/Python_Scripts/Objects/'

# Vocabularies
vocabulary_items = utils.load_yml_file(f'{objects_folder}/all_vocabularies.yml', 'r')
utils.create_new_vocabulary(session, vocabulary_items)

# Property types
property_items = utils.load_yml_file(f'{objects_folder}/all_property_types.yml', 'r')
utils.create_new_property_types(session, property_items)

# Objects
object_items = utils.load_yml_file(f'{objects_folder}/all_objects.yml', 'r')
utils.create_new_object_types(session, object_items)
