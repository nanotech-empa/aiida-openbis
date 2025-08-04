from langchain_core.tools import tool
from ai_agent import openbis_utils
from pydantic import BaseModel, Field, model_validator, field_validator
from typing import List, Dict, Literal, Annotated, Optional
from datetime import datetime
from schema.openbis_objects import Substance, Crystal, TwoDLayerMaterial

# Pydantic models for input validation
class TimestampInterval(BaseModel):
    begin_date_str: str = Field(
        ..., 
        description="Datetime in format YYYY-MM-DD HH:MM:SS, e.g., '2025-07-25 08:30:00'"
    )
    end_date_str: str = Field(
        ..., 
        description="Datetime in format YYYY-MM-DD HH:MM:SS, e.g., '2025-07-25 18:30:00'"
    )
    
    @model_validator(mode="after")
    def validate_datetime_order(self) -> "TimestampInterval":
        try:
            begin = datetime.strptime(self.begin_date_str, "%Y-%m-%d %H:%M:%S")
            end = datetime.strptime(self.end_date_str, "%Y-%m-%d %H:%M:%S")
        except ValueError as e:
            raise ValueError(f"Invalid datetime format: {e}")

        if begin >= end:
            raise ValueError("begin_date_str must be earlier than end_date_str")

        return self

# Generic openBIS tools

@tool("get_openbis_objects")
def get_openbis_objects(type: str, collection_identifier: str = None) -> List[Dict]:
    """
    Return list of openBIS objects based on the type.
    
    Args:
        type (str): openBIS object type
        collection_identifier (str): openBIS collection identifier
    """
    objects = openbis_utils.get_openbis_objects(type, collection_identifier)
    objects_data = []
    for obj in objects:
        obj_data = openbis_utils.get_openbis_object_data(obj)
        objects_data.append(obj_data)
    
    return objects_data[0:50] # TODO: Remove this restriction or adapt the function to ask the user if one really wants that

@tool("get_openbis_object_by_permId")
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

@tool("get_openbis_object_by_name")
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

@tool("get_openbis_object_by_date")
def get_openbis_objects_by_date(time_interval: TimestampInterval) -> List[Dict]:
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

@tool("get_openbis_object_by_description")
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

# Specific tools

# @tool("get_substances_by_attributes")
def get_substances_by_attributes(substance: Substance) -> List[Dict]:
    """
    Search for substances in openBIS using one or more identifying attributes.

    This tool retrieves substances by matching the given Empa identifier (empa number and batch together), e.g. 
    704a (empa number is 704 and batch is a), or any molecular descriptor (SMILES, sum formula, IUPAC name). 
    You can specify any subset of these fields. Unprovided fields will be ignored in the search.

    Args:
        substance (Substance): A data object describing the substance to search for. It may include:
            - empa_number (int, optional): Empa identifier number, e.g. 704.
            - batch (str, optional): Batch letter of the substance, e.g. 'a'.
            - molecules (List[Molecule], optional): A list of one or more molecules that describe the substance. 
              Each molecule may include:
                - smiles (str, optional): SMILES string, e.g. "CCO"
                - sum_formula (str, optional): Molecular sum formula, e.g. "CH4"
                - iupac_name (str, optional): IUPAC name, e.g. "benzene"

    Returns:
        List[Dict]: A list of dictionaries, where each dictionary represents a substance object from openBIS 
        matching the provided criteria.

    Example:
        >>> get_substances_by_attributes(
        ...     Substance(
        ...         empa_number=704,
        ...         batch="a",
        ...         molecules=[Molecule(smiles="CCO")]
        ...     )
        ... )
        # Returns all matching substance objects with Empa ID 704a and containing a molecule with SMILES "CCO".
    """
    
    obj_type = "SUBSTANCE"
    obj_type_lower = obj_type.lower()
    objects = openbis_utils.get_openbis_objects(obj_type)
    objects_data = []
    
    for obj in objects:
        obj_props = obj.props.all()
        obj_empa_number = obj_props[f"{obj_type_lower}.empa_number"]
        if obj_empa_number:
            obj_empa_number = int(obj_empa_number)
        else:
            obj_empa_number = 0
            
        obj_batch = obj_props[f"{obj_type_lower}.batch"]
        if obj_batch is None:
            obj_batch = ""
            
        obj_molecules = obj_props[f"{obj_type_lower}.molecules"]
        
        load_obj_data = True
        if substance.empa_number > 0:
            load_obj_data = obj_empa_number == substance.empa_number
        
        if substance.batch:
            load_obj_data = (obj_batch == substance.batch and load_obj_data)
        
        if substance.molecules:
            if obj_molecules:
                for molecule_permId in obj_molecules:
                    molecule_obj = openbis_utils.get_openbis_object(molecule_permId)
                    molecule_obj_props = molecule_obj.props.all()
                    molecule_smiles = molecule_obj_props[f"{obj_type_lower}.smiles"]
                    molecule_formula = molecule_obj_props[f"{obj_type_lower}.sum_formula"]
                    molecule_iupac = molecule_obj_props[f"{obj_type_lower}.iupac_name"]
                    
                    for prompt_molecule in substance.molecules:
                        if prompt_molecule.smiles:
                            load_obj_data = (prompt_molecule.smiles == molecule_smiles and load_obj_data)
                        
                        if prompt_molecule.sum_formula:
                            load_obj_data = (prompt_molecule.sum_formula == molecule_formula and load_obj_data)
                        
                        if prompt_molecule.iupac_name:
                            load_obj_data = (prompt_molecule.iupac_name == molecule_iupac and load_obj_data)
        
        if load_obj_data:
            obj_data = openbis_utils.get_openbis_object_data(obj)
            objects_data.append(obj_data)
                                
    return objects_data

@tool("get_crystals_by_attributes")
def get_crystals_by_attributes(crystal: Crystal) -> List[Dict]:
    """
    Search for crystals in openBIS using one or more identifying attributes.

    This tool retrieves crystals by matching the given crystal concept descriptor (face, material). 
    You can specify any subset of these fields. Unprovided fields will be ignored in the search.

    Args:
        crystal (Crystal): A data object describing the crystal to search for. It may include:
            - crystal_concept (CrystalConcept, optional): crystal concept that contains properties of the crystal
              Each crystal_concept may include:
                - face (str, optional): Material face, e.g. 111 or 10000
                - material (str, optional): Crystal material, e.g. Au, Pt, Ag, Pd

    Returns:
        List[Dict]: A list of dictionaries, where each dictionary represents a crystal object from openBIS 
        matching the provided criteria.

    Example:
        >>> get_crystals_by_attributes(
        ...     Crystal(
        ...         crystal_concepts=[CrystalConcept(face="111", material="Au")]
        ...     )
        ... )
        # Returns all matching crystal objects containing a crystal concept with face 111 and material Au.
        
        If the user asks for crystals Au111 or Gold-111, its means that one wants crystals with a crystal 
        concept with material Au and face 111
    """
    obj_type = "CRYSTAL"
    obj_type_lower = obj_type.lower()
    objects = openbis_utils.get_openbis_objects(obj_type)
    objects_data = []
    
    for obj in objects:
        obj_props = obj.props.all()
        crystal_concept_permId = obj_props[f"{obj_type_lower}.concept"]
        load_obj_data = True
        if crystal.crystal_concept:
            if crystal_concept_permId:
                crystal_concept_obj = openbis_utils.get_openbis_object(crystal_concept_permId)
                crystal_concept_obj_props = crystal_concept_obj.props.all()
                crystal_concept_face = crystal_concept_obj_props[f"{obj_type_lower}.face"]
                crystal_concept_material = crystal_concept_obj_props[f"{obj_type_lower}.material"]
                
                if crystal.crystal_concept.face:
                    load_obj_data = (crystal.crystal_concept.face == crystal_concept_face and load_obj_data)
                
                if crystal.crystal_concept.material:
                    load_obj_data = (crystal.crystal_concept.material == crystal_concept_material and load_obj_data)
        
        if load_obj_data:
            obj_data = openbis_utils.get_openbis_object_data(obj)
            objects_data.append(obj_data)
                                
    return objects_data

@tool("get_2d_materials_by_attributes")
def get_2d_materials_by_attributes(two_d_material: TwoDLayerMaterial) -> List[Dict]:
    """
    Search for 2D layer materials in openBIS using one or more identifying attributes.

    This tool retrieves 2D layer materials by matching the given top layer material, layer count,
    substrate, growth/fabrication method.
    
    You can specify any subset of these fields. Unprovided fields will be ignored in the search.

    Args:
        two_d_material (TwoDLayerMaterial): A data object describing the 2d layer material 
        to search for. It may include:
            - top_layer_material (str, optional): Top layer material
            - layer_count (str, optional): Number of layers in string
            - substrate (str, optional): Substrate material
            - growth_method (str, optional): Growth/Fabrication method

    Returns:
        List[Dict]: A list of dictionaries, where each dictionary represents a 2d layer material object 
        from openBIS matching the provided criteria.

    Example:
        >>> get_2d_materials_by_attributes(
        ...     TwoDLayerMaterial(
        ...         top_layer_material='MoS2'
        ...     )
        ... )
        # Returns all matching 2d later material objects
    """
    
    obj_type = "2D_LAYER_MATERIAL"
    obj_type_lower = obj_type.lower()
    objects = openbis_utils.get_openbis_objects(obj_type)
    objects_data = []
    
    for obj in objects:
        obj_props = obj.props.all()
        obj_top_layer_material = obj_props[f"{obj_type_lower}.top_layer_material"]
        obj_layer_count = obj_props[f"{obj_type_lower}.layer_count"]
        obj_substrate = obj_props[f"{obj_type_lower}.substrate"]
        obj_growth_method = obj_props[f"{obj_type_lower}.growth_method"]
        
        load_obj_data = True
        if two_d_material.top_layer_material:
            load_obj_data = obj_top_layer_material == two_d_material.top_layer_material
        
        if two_d_material.layer_count:
            load_obj_data = (obj_layer_count == two_d_material.layer_count and load_obj_data)
            
        if two_d_material.substrate:
            load_obj_data = (obj_substrate == two_d_material.substrate and load_obj_data)
        
        if two_d_material.growth_method:
            load_obj_data = (obj_growth_method == two_d_material.growth_method and load_obj_data)
        
        if load_obj_data:
            obj_data = openbis_utils.get_openbis_object_data(obj)
            objects_data.append(obj_data)
                                
    return objects_data