# Load libraries
from aiida_openbis.utils import bisutils
import utils

# Connect to openBIS
session = bisutils.log_in(bisurl="openbis", bisuser="admin", bispasswd="changeit")

# Load property types and object types from YAML file
filepath = './publication_object.yml'
mode = 'r'
items = utils.load_yml_file(filepath, mode)

# Create new vocabulary
if 'vocabulary' in items:
    utils.create_new_vocabulary(session, items)

# Create new property types
if 'property_types' in items:
    utils.create_new_property_types(session, items)

# Create new object types
if 'object_types' in items:
    utils.create_new_object_types(session, items)