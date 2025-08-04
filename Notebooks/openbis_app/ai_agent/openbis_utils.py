import src.utils as utils
from functools import lru_cache

CONFIG_ELN = utils.get_aiidalab_eln_config()
OPENBIS_SESSION, SESSION_DATA = utils.connect_openbis(CONFIG_ELN["url"], CONFIG_ELN["token"])

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

@lru_cache(maxsize=100)
def get_openbis_objects(type = None, collection = None):
    if type:
        if collection:
            return OPENBIS_SESSION.get_objects(type = type, collection = collection)
        else:
            return OPENBIS_SESSION.get_objects(type = type)
    else:
        return OPENBIS_SESSION.get_objects()

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
    if obj_parents:
        obj_parents = obj_parents.df.permId.to_list()
        
    obj_parents_data = []
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