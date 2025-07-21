import src.utils as utils
import asyncio
from zoneinfo import ZoneInfo
from datetime import datetime
from functools import lru_cache
from langchain_core.tools import tool
from langgraph.checkpoint.memory import InMemorySaver
from langgraph.prebuilt import create_react_agent
from langchain_google_genai import ChatGoogleGenerativeAI
from typing import List, Dict

# OpenBIS configuration and session setup
CONFIG_ELN = utils.get_aiidalab_eln_config()
OPENBIS_SESSION, SESSION_DATA = utils.connect_openbis(CONFIG_ELN["url"], CONFIG_ELN["token"])

def get_current_time() -> str:
    """
    Get the current time formatted nicely for display.
    """
    # Get current datetime
    now = datetime.now(ZoneInfo("Europe/Zurich"))
    # Format as a human-readable sentence
    formatted = now.strftime("%A, %B %d, %Y at %H:%M")
    return f"It is currently {formatted}."

@lru_cache(maxsize=5000)
def _get_openbis_object(permId):
    try:
        obj = OPENBIS_SESSION.get_object(permId)
        return obj
    except ValueError as e:
        return None

@lru_cache(maxsize=5000)
def _get_openbis_property_type(key):
    try:
        prop = OPENBIS_SESSION.get_property_type(key)
        return prop
    except ValueError as e:
        return None

@lru_cache(maxsize=100)
def _get_openbis_objects_by_type(type):
    openbis_objs = OPENBIS_SESSION.get_objects(type = type)
    objs = {}
    for openbis_obj in openbis_objs:
        objs[openbis_obj.permId] = openbis_obj
    return objs

@lru_cache(maxsize=1000)
def _get_openbis_experiments():
    experiments = {}
    openbis_experiments = OPENBIS_SESSION.get_experiments(type = "EXPERIMENT")
    for exp in openbis_experiments:
        experiments[exp.permId] = exp
    return experiments

@lru_cache(maxsize=1000)
def _get_openbis_experiment(permId):
    return OPENBIS_SESSION.get_experiment(permId)

@lru_cache(maxsize=1000)
def _get_substance_data(permId):
    substance_data = {}
    obj = _get_openbis_object(permId)
    if obj:
        substance_properties = obj.props.all()
        for key, value in substance_properties.items():
            if value:
                prop_type = _get_openbis_property_type(key)
                prop_datatype = prop_type.dataType
                if prop_datatype == "SAMPLE":
                    if key == "molecules":
                        molecules_data = []
                        for molecule in value:
                            molecule = _get_openbis_object(molecule)
                            if molecule:
                                molecule_data = {}
                                molecule_properties = molecule.props.all()
                                for molecule_key, molecule_value in molecule_properties.items():
                                    if molecule_value:
                                        molecule_data[molecule_key] = molecule_value
                                
                                molecules_data.append(molecule_data)
                        
                        substance_data["molecules"] = molecules_data
                    else:
                        if isinstance(value, list):
                            prop_data = []
                            for obj_permId in value:
                                openbis_object = _get_openbis_object(obj_permId)
                                if openbis_object:
                                    openbis_object_properties = openbis_object.props.all()
                                    openbis_name = openbis_object_properties.get("name", "Unknown")
                                    prop_data.append({"name": openbis_name})
                            
                            substance_data[key] = prop_data
                
                else:
                    substance_data[key] = value
    
    return substance_data

@lru_cache(maxsize=100)
def _get_experiment_data(permId):
    experiment = _get_openbis_experiment(permId)
    exp_objects = experiment.get_objects()
    experiment_data = {
        "description": "",
        "date": "",
        "name": experiment.props["name"],
        "samples": [],
        "substances": [],
        "crystals": [],
        "instruments": []
    }
    experiment_info = ""
    for obj in exp_objects:
        obj = _get_openbis_object(obj.permId)
        if obj.type == "PREPARATION":
            preparation_properties = obj.props.all()
            preparation_info = "Preparation: "
            for prep_key, prep_value in preparation_properties.items():
                if prep_value:
                    prop_type = _get_openbis_property_type(prep_key)
                    prop_label = prop_type.label
                    preparation_info += f"{prop_label}: {prep_value}\n"
            preparation_children = obj.children
            
            for child in preparation_children:
                child_obj = _get_openbis_object(child)
                if child_obj.type == "PROCESS_STEP":
                    step_properties = child_obj.props.all()
                    step_info = "Process step: "
                    
                    for step_key, step_value in step_properties.items():
                        if step_value:
                            if step_key == "actions":
                                actions_objs = []
                                actions_info = "Actions: "
                                for action_permId in step_value:
                                    action_obj = _get_openbis_object(action_permId)
                                    action_type = action_obj.type
                                    action_properties = action_obj.props.all()
                                    actions_info += f"{action_type}: "
                                    for action_key, action_value in action_properties.items():
                                        if action_value:
                                            if action_key == "substance":
                                                substance_info = "Substance: "
                                                substance_obj = _get_openbis_object(action_value)
                                                substance_props = substance_obj.props.all()
                                                prop_type = _get_openbis_property_type(action_key)
                                                prop_label = prop_type.label
                                                prop_datatype = prop_type.dataType
                                                if prop_datatype == "SAMPLE":
                                                    for substance_key, substance_value in substance_props.items():
                                                        if substance_key == "molecules":
                                                            molecule_info = "Molecules: "
                                                            for molecule in substance_value:
                                                                molecule = _get_openbis_object(molecule)
                                                                if molecule:
                                                                    molecule_properties = molecule.props.all()
                                                                    for molecule_key, molecule_value in molecule_properties.items():
                                                                        if molecule_value:
                                                                            molecule_prop_type = _get_openbis_property_type(molecule_key)
                                                                            molecule_prop_label = molecule_prop_type.label
                                                                            molecule_info += f"{molecule_prop_label}: {molecule_value}\n"
                                                            
                                                            substance_info += molecule_info
                                                            
                                                        else:
                                                            prop_type = _get_openbis_property_type(substance_key)
                                                            prop_label = prop_type.label
                                                            substance_info += f"{prop_label}: {substance_value}\n"
                                                
                                                actions_info += f"{substance_info}"
                                            
                                            elif action_key == "component":
                                                component_obj = _get_openbis_object(action_value)
                                                component_properties = component_obj.props.all()
                                                component_name = component_properties.get("name", "Unknown")
                                                step_info += f"Component: {component_name}\n"
                                                
                                            else:
                                                prop_type = _get_openbis_property_type(action_key)
                                                prop_label = prop_type.label
                                                actions_info += f"{prop_label}: {action_value}\n"
                                step_info += actions_info
                                    
                            elif step_key == "observables":
                                pass
                            elif step_key == "instrument":
                                instrument_obj = _get_openbis_object(step_value)
                                instrument_properties = instrument_obj.props.all()
                                instrument_name = instrument_properties.get("name", "Unknown")
                                step_info += f"Instrument: {instrument_name}\n"
                            else:
                                prop_type = _get_openbis_property_type(step_key)
                                prop_label = prop_type.label
                                step_info += f"{prop_label}: {step_value}\n"
            
                preparation_info += step_info
            
            experiment_info += preparation_info
            
        elif obj.type == "MEASUREMENT_SESSION":
            pass
    
    all_experiments_info += f"Experiment: {experiment_info}"
    
    return all_experiments_info

# Converts a Python function into a LangChain-compatible tool 
# that the agent can call automatically
@tool("get_substances")
def _get_substances() -> List[Dict]:
    """Returns information about all substances or molecules available in the inventory of openBIS."""
    substances_data = []
    substances_objects = _get_openbis_objects_by_type("SUBSTANCE")
    for obj_permId, obj in substances_objects.items():
        substance_data = _get_substance_data(obj_permId)
        substances_data.append(substance_data)

    return substances_data

@tool("get_substances_by_name")
def _get_substances_by_name(name: str) -> List[Dict]:
    """Returns substance or molecule metadata available in the inventory of openBIS by name."""
    substances_objects = _get_openbis_objects_by_type("SUBSTANCE")
    substances_data = []
    for obj_permId, obj in substances_objects.items():
        substance_properties = obj.props.all()
        if substance_properties["name"] == name:
            substance_data = _get_substance_data(obj_permId)
            substances_data.append(substance_data)

    return substances_data
        
@tool("get_substances_by_smiles")
def _get_substances_by_smiles(smiles: str) -> List[Dict]:
    """Returns substance or molecule metadata available in the inventory of openBIS by SMILES string."""
    substances_data = []
    substances_objects = _get_openbis_objects_by_type("SUBSTANCE")
    for substance_obj_permId, substance_obj in substances_objects.items():
        smiles_found = False
        substance_properties = substance_obj.props.all()
        for substance_prop_key, substance_prop_value in substance_properties.items():
            if substance_prop_key == "molecules" and substance_prop_value:
                for molecule_permId in substance_prop_value:
                    molecule_obj = _get_openbis_object(molecule_permId)
                    if molecule_obj:
                        molecule_properties = molecule_obj.props.all()
                        if molecule_properties["smiles"] == smiles:
                            smiles_found = True
                
                if smiles_found:
                    substance_data = _get_substance_data(substance_obj_permId)
                    substances_data.append(substance_data)

    return substances_data

@tool("get_substances_by_empa_identifier")
def _get_substances_by_empa_identifier(empa_number: str, batch: str) -> List[Dict]:
    """Returns substance or molecule metadata available in the inventory of openBIS by Empa number and batch, e.g., molecule 125a or 700b."""
    substances_objects = _get_openbis_objects_by_type("SUBSTANCE")
    substances_data = []
    for obj_permId, obj in substances_objects.items():
        substance_properties = obj.props.all()
        if substance_properties["empa_number"] == empa_number and substance_properties["batch"] == batch:
            substance_data = _get_substance_data(obj_permId)
            substances_data.append(substance_data)
            break

    return substances_data

@tool("get_crystals")
def _get_crystals() -> str:
    """Returns information about all crystals available in the inventory of openBIS."""
    crystal_objects = _get_openbis_objects_by_type("CRYSTAL")
    crystals_info = ""
    for obj_id, obj in crystal_objects.items():
        crystal_properties = obj.props.all()
        for key, value in crystal_properties.items():
            if value:
                prop_type = _get_openbis_property_type(key)
                prop_label = prop_type.label
                prop_datatype = prop_type.dataType
                if prop_datatype == "SAMPLE":
                    if key == "crystal_concept":
                        crystal_concept = _get_openbis_object(value)
                        if crystal_concept:
                            crystal_concept_properties = crystal_concept.props.all()
                            value = ""
                            for concept_key, concept_value in crystal_concept_properties.items():
                                if concept_value:
                                    concept_prop_type = _get_openbis_property_type(key)
                                    concept_prop_label = concept_prop_type.label
                                    value += f"{concept_prop_label}: {concept_value}\n"
                    else:
                        openbis_object = _get_openbis_object(value)
                        if openbis_object:
                            openbis_object_properties = openbis_object.props.all()
                            openbis_name = openbis_object_properties.get("name", "Unknown")
                            value = f"Name: {openbis_name}\n"

                crystals_info += f"{prop_label}: {value}\n"

    return crystals_info

@tool("get_experiments")
def _get_experiments() -> List[Dict]:
    """Returns a list of experiments, each with all the preparations, process steps, actions, observables and materials used.
    The materials include crystals, 2D materials, substances (molecules), among others."""
    
    experiments = _get_openbis_experiments()
    experiments_data = []
    for exp_permId, exp in experiments.items():
        experiment_data = {
            "preparations": [],
            "measurements": []
        }
        exp_objects = exp.get_objects()
        for obj in exp_objects:
            obj = _get_openbis_object(obj.permId)
            if obj.type == "PREPARATION":
                preparation_data = {
                    "name": "",
                    "process_steps": []
                }
                preparation_properties = obj.props.all()
                preparation_data["name"] = preparation_properties["name"]
                preparation_children = obj.children
                
                for child in preparation_children:
                    child_obj = _get_openbis_object(child)
                    if child_obj.type == "PROCESS_STEP":
                        process_step_properties = child_obj.props.all()
                        process_step_data = {
                            "name": "",
                            "actions": [],
                            "observables": [],
                            "instrument": "",
                            "sample_in": "",
                            "sample_out": ""
                        }
                        process_step_data["name"] = process_step_properties["name"]
                        
                        # Get sample in name
                        process_step_parents = child_obj.parents
                        for parent in process_step_parents:
                            parent_obj = _get_openbis_object(parent)
                            if parent_obj.type == "SAMPLE":
                                process_step_data["sample_in"] = parent_obj.props["name"]
                                break
                        
                        # Get sample out name
                        process_step_children = child_obj.children
                        for process_step_child in process_step_children:
                            process_step_child_obj = _get_openbis_object(process_step_child)
                            if process_step_child_obj.type == "SAMPLE":
                                process_step_data["sample_out"] = process_step_child_obj.props["name"]
                                break
                        
                        for step_key, step_value in process_step_properties.items():
                            if step_value:
                                if step_key == "actions":
                                    for action_permId in step_value:
                                        action_data = {
                                            "name": "",
                                            "type": "",
                                            "description": "",
                                            "duration": "",
                                            "target_temperature": "",
                                            "substrate_temperature": "",
                                            "substance": "",
                                            "component": "",
                                            "component_settings": "",
                                            "dosing_gas": "",
                                            "pressure": "",
                                            "sputter_ion": "",
                                            "current": "",
                                            "angle": "",
                                            "comments": ""
                                        }
                                        
                                        action_obj = _get_openbis_object(action_permId)
                                        action_properties = action_obj.props.all()
                                        action_data["name"] = action_properties.pop("name")
                                        action_data["type"] = action_obj.type.code
                                        
                                        for action_key, action_value in action_properties.items():
                                            if action_value:
                                                prop_type = _get_openbis_property_type(action_key)
                                                prop_datatype = prop_type.dataType
                                                if prop_datatype == "SAMPLE":
                                                    prop_sample = _get_openbis_object(action_value)
                                                    prop_sample_name = prop_sample.props["name"]
                                                    action_data[action_key] = prop_sample_name
                                                else:
                                                    action_data[action_key] = action_value
                                                    
                                        process_step_data["actions"].append(action_data)
                                        
                                elif step_key == "observables":
                                    for observable_permId in step_value:
                                        observable_data = {
                                            "name": "",
                                            "type": "",
                                            "description": "",
                                            "duration": "",
                                            "channel_name": "",
                                            "component": "",
                                            "component_settings": "",
                                            "comments": ""
                                        }
                                        
                                        observable_obj = _get_openbis_object(observable_permId)
                                        observable_properties = observable_obj.props.all()
                                        observable_data["name"] = observable_properties.pop("name")
                                        observable_data["type"] = observable_obj.type.code
                                        
                                        for observable_key, observable_value in observable_properties.items():
                                            if observable_value:
                                                prop_type = _get_openbis_property_type(observable_key)
                                                prop_datatype = prop_type.dataType
                                                if prop_datatype == "SAMPLE":
                                                    prop_sample = _get_openbis_object(observable_value)
                                                    prop_sample_name = prop_sample.props["name"]
                                                    observable_data[observable_key] = prop_sample_name
                                                else:
                                                    observable_data[observable_key] = observable_value
                                                    
                                        process_step_data["observables"].append(observable_data)
                                elif step_key == "instrument":
                                    instrument_obj = _get_openbis_object(step_value)
                                    process_step_data["instrument"] = instrument_obj.props["name"]
                                else:
                                    process_step_data[step_key] = step_value

                        preparation_data["process_steps"].append(process_step_data)
                
                experiment_data["preparations"].append(preparation_data)
                
            elif obj.type == "MEASUREMENT_SESSION":
                pass
        
        experiments_data.append(experiment_data)
    
    return experiments_data

class OpenBISAgent():
    """
    An agent that can answer questions about openBIS inventories.
    It uses a LangChain model and tools to interact with the openBIS session.
    """
    def __init__(self, google_api_key):
        self.google_api_key = google_api_key
        self.llm_model = ChatGoogleGenerativeAI(model = "models/gemini-2.5-flash",google_api_key = self.google_api_key)
        self.system_prompt = """
            You are an helpful assistant that can answer questions about experiments, simulations, and all the inventories inside 
            openBIS. {get_current_time()}
        """
        
        self._tools = [
            _get_substances,
            _get_substances_by_name,
            _get_substances_by_smiles,
            _get_substances_by_empa_identifier,
            _get_crystals,
            _get_experiments
        ]
        
        self.ai_agent = create_react_agent(
            model = self.llm_model,
            tools = self._tools,
        )

    def ask_question(self, user_prompt: str):
        """
        Ask a question to the agent and get a response.
        
        Args:
            user_prompt (str): The question to ask the agent.
        
        Returns:
            str: The response from the agent.
        """
        # Builds a LangGraph agent that interleaves reasoning and tool execution (ReAct pattern)
        # Runs the full workflow, returning a state object 
        # where the last message contains the model’s response after any tool invocations
        response = self.ai_agent.invoke({
            "messages": [
                {"role": "system", "content": self.system_prompt},
                {"role": "user", "content": user_prompt}
            ]
        })

        return response["messages"]