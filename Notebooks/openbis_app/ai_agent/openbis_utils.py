import sys
import os

current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)
sys.path.append(parent_dir)

from src import utils
from functools import lru_cache

OPENBIS_SESSION, SESSION_DATA = utils.connect_openbis_aiida()

@lru_cache(maxsize=5000)
def get_openbis_property_type(key):
    try:
        return OPENBIS_SESSION.get_property_type(key)
    except ValueError as e:
        return None

@lru_cache(maxsize=5000)
def get_openbis_object(permId):
    try:
        return OPENBIS_SESSION.get_object(permId)
    except ValueError as e:
        return None

def get_openbis_objects(**kwargs):
    return OPENBIS_SESSION.get_objects(**kwargs)

@lru_cache(maxsize=5000)
def get_openbis_object_data(obj, depth=1):
    obj_props = obj.props.all()
        
    complete_obj_props = {}
    for prop_key, prop_value in obj_props.items():
        if prop_value:
            prop_type = get_openbis_property_type(prop_key)
            if prop_type.dataType == "SAMPLE" and depth > 0:
                if isinstance(prop_value, list):
                    complete_prop_value = []
                    for prop_obj_permId in prop_value:
                        prop_obj = get_openbis_object(prop_obj_permId)
                        prop_obj_data = get_openbis_object_data(prop_obj, depth = depth - 1)
                        complete_prop_value.append(prop_obj_data)
                        
                    prop_value = complete_prop_value
                else:
                    prop_obj = get_openbis_object(prop_value)
                    prop_obj_data = get_openbis_object_data(prop_obj, depth = depth - 1)
                    prop_value = prop_obj_data
        
        complete_obj_props[prop_key] = prop_value
    
    obj_parents = obj.get_parents()
    obj_parents_data = []
    if obj_parents:
        obj_parents = obj_parents.df.permId.to_list()
        
        for parent in obj_parents:
            parent_obj = get_openbis_object(parent)
            parent_obj_data = get_openbis_object_data(parent_obj, depth = depth - 1)
            obj_parents_data.append(parent_obj_data)
        
    obj_data = {
        "permId": obj.permId,
        "type": obj.type.code,
        "properties": complete_obj_props,
        "parents": obj_parents_data,
        "registration_date": obj.registrationDate
    }
    
    return obj_data

@lru_cache(maxsize=1000)
def get_openbis_experiments():
    return OPENBIS_SESSION.get_experiments(type = "EXPERIMENT")

@lru_cache(maxsize=1000)
def get_openbis_experiment(permId):
    try:
        return OPENBIS_SESSION.get_experiment(permId)
    except ValueError as e:
        return None