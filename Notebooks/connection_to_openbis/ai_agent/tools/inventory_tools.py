from langchain_core.tools import tool
from ai_agent import openbis_utils
from pydantic import BaseModel, Field, model_validator, field_validator
from typing import List, Dict, Literal, Annotated, Optional
from datetime import datetime

@tool("get_openbis_objects")
def get_openbis_objects(type: str, collection_identifier: str = None):
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
    
    return objects_data

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

class SubstanceEmpaIdentifier(BaseModel):
    empa_number: int = Field(
        ..., 
        description = "Integer value given to new substances inside Empa"
    )
    batch: str = Field(
        ..., 
        description = "Letter given to the batch of substances, e.g., a for the first, b for the second."
    )

@tool("get_substance_by_empa_id")
def get_substances_by_empa_id(empa_identifier: SubstanceEmpaIdentifier) -> List[Dict]:
    """
    Return a list of substances available in openBIS by empa number and batch.
    Example:
        - User: Is substance 704a available? empa_number is 704 and batch is a.
    
    Args:
        empa_identifier: An object containing:
            - empa_number (int): Empa identifier given to substances in the labs
            - batch (str): Letter given to the different batches of a substance in the labs
    
    Returns:
        A list of dictionaries with openBIS objects matching the Empa identifier.
    """
    obj_type = "SUBSTANCE"
    objects = openbis_utils.get_openbis_objects(obj_type)
    objects_data = []
    for obj in objects:
        obj_empa_number = obj.props["empa_number"]
        if obj_empa_number:
            obj_empa_number = int(obj_empa_number)
            
        if obj_empa_number == empa_identifier.empa_number and obj.props["batch"] == empa_identifier.batch:
            obj_data = openbis_utils.get_openbis_object_data(obj)
            objects_data.append(obj_data)
    
    return objects_data

@tool("get_substance_by_smiles")
def get_substances_by_smiles(smiles: str) -> List[Dict]:
    """
    Return a list of substances, molecules, or chemicals available in openBIS by SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        A list of dictionaries with openBIS objects matching the SMILES string, e.g. CCO.
    """
    obj_type = "SUBSTANCE"
    objects = openbis_utils.get_openbis_objects(obj_type)
    objects_data = []
    for obj in objects:
        obj_molecules = obj.props["molecules"]
        if obj_molecules:
            for molecule_permId in obj_molecules:
                molecule_obj = openbis_utils.get_openbis_object(molecule_permId)
                molecule_smiles = molecule_obj.props["smiles"]
                if molecule_smiles:
                    if molecule_smiles.lower() == smiles.lower():
                        obj_data = openbis_utils.get_openbis_object_data(obj)
                        objects_data.append(obj_data)
    
    return objects_data

@tool("get_substance_by_formula")
def get_substances_by_formula(chemical_formula: str) -> List[Dict]:
    """
    Return a list of substances, molecules, or chemicals available in openBIS by chemical formula.
    
    Args:
        chemical_formula (str): Chemical formula of the molecule, e.g. CH4
    
    Returns:
        A list of dictionaries with openBIS objects matching the chemical formula.
    """
    obj_type = "SUBSTANCE"
    objects = openbis_utils.get_openbis_objects(obj_type)
    objects_data = []
    for obj in objects:
        obj_molecules = obj.props["molecules"]
        if obj_molecules:
            for molecule_permId in obj_molecules:
                molecule_obj = openbis_utils.get_openbis_object(molecule_permId)
                molecule_formula = molecule_obj.props["sum_formula"]
                if molecule_formula:
                    if molecule_formula.lower() == chemical_formula.lower():
                        obj_data = openbis_utils.get_openbis_object_data(obj)
                        objects_data.append(obj_data)
    
    return objects_data

class MoleculeObject(BaseModel):
    smiles: Optional[str] = Field(default="", description = "SMILES string for the substance, e.g. CCO")
    sum_formula: Optional[str] = Field(default="", description = "Molecular sum formula, e.g. CH4")
    iupac_name: Optional[str] = Field(default="", description="IUPAC standardised chemical name")

class SubstanceObject(BaseModel):
    empa_number: Optional[int] = Field(default=0, description = "Integer value given to new substances inside Empa")
    batch: Optional[str] = Field(default="", description = "Letter given to the batch of substances, e.g., a for the first, b for the second.")
    molecules: Optional[List[MoleculeObject]] = Field(
        default_factory=list, 
        description = "List of molecules that substance contains",
    )
    
    class Config:
        # This ensures proper JSON schema generation
        json_schema_extra = {
            "properties": {
                "molecules": {
                    "type": "array",
                    "items": {
                        "type": "object",
                        "properties": {
                            "smiles": {"type": "string", "description": "SMILES string for the substance, e.g. CCO"},
                            "sum_formula": {"type": "string", "description": "Molecular sum formula, e.g. CH4"},
                            "iupac_name": {"type": "string", "description": "IUPAC standardised chemical name"}
                        }
                    }
                }
            }
        }

# @tool("get_substances_by_attributes")
def get_substances_by_attributes(substance: SubstanceObject) -> List[Dict]:
    """
    Search for substances in openBIS using one or more identifying attributes.

    This tool retrieves substances by matching the given Empa identifier or any molecular descriptor 
    (SMILES, sum formula, IUPAC name). You can specify any subset of these fields—unprovided fields 
    will be ignored in the search.

    Args:
        substance (SubstanceObject): A data object describing the substance to search for. It may include:
            - empa_number (int, optional): Empa identifier number, e.g. 704.
            - batch (str, optional): Batch letter of the substance, e.g. 'a'.
            - molecules (List[MoleculeObject], optional): A list of one or more molecules that describe the substance. 
              Each molecule may include:
                - smiles (str, optional): SMILES string, e.g. "CCO"
                - sum_formula (str, optional): Molecular sum formula, e.g. "CH4"
                - iupac_name (str, optional): IUPAC name, e.g. "benzene"

    Returns:
        List[Dict]: A list of dictionaries, where each dictionary represents a substance object from openBIS 
        matching the provided criteria.

    Example:
        >>> get_substances_by_attributes(
        ...     SubstanceObject(
        ...         empa_number=704,
        ...         batch="a",
        ...         molecules=[MoleculeObject(smiles="CCO")]
        ...     )
        ... )
        # Returns all matching substance objects with Empa ID 704a and containing a molecule with SMILES "CCO".
    """
    # return []
    obj_type = "SUBSTANCE"
    objects = openbis_utils.get_openbis_objects(obj_type)
    objects_data = []
    
    for obj in objects:
        obj_props = obj.props.all()
        obj_empa_number = obj_props["empa_number"]
        if obj_empa_number:
            obj_empa_number = int(obj_empa_number)
        else:
            obj_empa_number = 0
            
        obj_batch = obj_props["batch"]
        if obj_batch is None:
            obj_batch = ""
            
        obj_molecules = obj_props["molecules"]
        
        if substance.empa_number:
            if substance.batch:
                if obj_empa_number == substance.empa_number and obj_batch == substance.batch:
                    obj_data = openbis_utils.get_openbis_object_data(obj)
                    if obj_data not in objects_data:
                        objects_data.append(obj_data)
            else:
                if obj_empa_number == substance.empa_number:
                    obj_data = openbis_utils.get_openbis_object_data(obj)
                    if obj_data not in objects_data:
                        objects_data.append(obj_data)
        
        if substance.molecules:
            if obj_molecules:
                for molecule_permId in obj_molecules:
                    molecule_obj = openbis_utils.get_openbis_object(molecule_permId)
                    molecule_obj_props = molecule_obj.props.all()
                    molecule_smiles = molecule_obj_props["smiles"]
                    molecule_formula = molecule_obj_props["sum_formula"]
                    molecule_iupac = molecule_obj_props["iupac_name"]
                    
                    for prompt_molecule in substance.molecules:
                        prompt_molecule_smiles = prompt_molecule.smiles
                        prompt_molecule_formula = prompt_molecule.sum_formula
                        prompt_molecule_iupac = prompt_molecule.iupac_name
                        
                        if prompt_molecule_smiles == molecule_smiles:
                            obj_data = openbis_utils.get_openbis_object_data(obj)
                            if obj_data not in objects_data:
                                objects_data.append(obj_data)
                        
                        if prompt_molecule_formula == molecule_formula:
                            obj_data = openbis_utils.get_openbis_object_data(obj)
                            if obj_data not in objects_data:
                                objects_data.append(obj_data)
                        
                        if prompt_molecule_iupac == molecule_iupac:
                            obj_data = openbis_utils.get_openbis_object_data(obj)
                            if obj_data not in objects_data:
                                objects_data.append(obj_data)
                                
    return objects_data