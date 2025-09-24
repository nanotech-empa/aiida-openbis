from langchain_core.tools import tool
from typing import List, Dict
from datetime import datetime

import sys
sys.path.append('../ai_agent')

import openbis_utils
import utils

@tool
def get_openbis_objects(obj_type: str, collection_identifier: str = None) -> List[Dict]:
    """
    Return list of openBIS objects based on the type.
    
    Args:
        obj_type (str): openBIS object type
        collection_identifier (str): openBIS collection identifier
    
    E.g.:
        User: Give me all instruments in openBIS
        Tool: get_openbis_objects(obj_type="INSTRUMENT")
        User: Give me all instruments in collection /INSTRUMENTS/STM
        Tool: get_openbis_objects(obj_type="INSTRUMENT", collection_identifier="/INSTRUMENTS/STM")
    """
    if collection_identifier:
        objects = openbis_utils.get_openbis_objects(type = obj_type, collection = collection_identifier)
    else:
        objects = openbis_utils.get_openbis_objects(type = obj_type)
    
    objects_data = []
    for obj in objects:
        obj_data = openbis_utils.get_openbis_object_data(obj)
        objects_data.append(obj_data)
    
    return objects_data[0:100] # TODO: Remove this restriction or adapt the function to ask the user if one really wants that

@tool
def get_openbis_object_by_permId(permId: str) -> Dict:
    """
    Return an object available in openBIS by permanent ID (permId). PermId follows the following 
    pattern: YYYYMMDDHHMMSSmmm-####. Example: 20250717115421436-984.
    
    Args:
        permId (str): openBIS object permanent ID (permId). Follows the pattern YYYYMMDDHHMMSSmmm-####.
        Example: 20250717115421436-984
    
    """
    obj = openbis_utils.get_openbis_object(permId)
    obj_data = openbis_utils.get_openbis_object_data(obj)
    
    return obj_data

@tool
def get_openbis_objects_by_name(name: str) -> List[Dict]:
    """
    Return objects available in openBIS by name.
    
    Args:
        name (str): openBIS object name
    """
    objects = openbis_utils.get_openbis_objects()
    objects_data = []
    for obj in objects:
        if obj.props["name"] == name:
            obj_data = openbis_utils.get_openbis_object_data(obj)
            objects_data.append(obj_data)
    
    return objects_data

@tool
def get_openbis_objects_by_date(time_interval: utils.TimestampInterval) -> List[Dict]:
    """
    Returns objects available in openBIS that were registered between the given time interval.

    Args:
        time_interval: An object containing:
            - begin_date_str (str): Start of the time interval in format YYYY-MM-DD HH:MM:SS
            - end_date_str (str): End of the time interval in format YYYY-MM-DD HH:MM:SS

    Returns:
        A list of dictionaries with openBIS objects matching the interval.
    """
    begin_date = datetime.strptime(time_interval.begin_date_str, "%Y-%m-%d %H:%M:%S")
    end_date = datetime.strptime(time_interval.end_date_str, "%Y-%m-%d %H:%M:%S")
    
    objects = openbis_utils.get_openbis_objects()
    objects_data = []
    for obj in objects:
        obj_registration_date = datetime.strptime(obj.registrationDate, "%Y-%m-%d %H:%M:%S")
        if begin_date <= obj_registration_date and end_date >= obj_registration_date:
            obj_data = openbis_utils.get_openbis_object_data(obj)
            objects_data.append(obj_data)
    
    return objects_data

@tool
def get_openbis_objects_by_description(description: str) -> List[Dict]:
    """
    Return objects available in openBIS by description.
    """
    objects = openbis_utils.get_openbis_objects()
    objects_data = []
    for obj in objects:
        obj_props = obj.props.all()
        obj_description = obj_props.get("description", "")
        # TODO: Here the tool needs to check by similarity and not by equality
        if obj_description == description:
            obj_data = openbis_utils.get_openbis_object_data(obj)
            objects_data.append(obj_data)
    
    return objects_data
