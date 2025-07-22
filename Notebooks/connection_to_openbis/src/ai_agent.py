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
from inspect import getmembers

# OpenBIS configuration and session setup
CONFIG_ELN = utils.get_aiidalab_eln_config()
OPENBIS_SESSION, SESSION_DATA = utils.connect_openbis(CONFIG_ELN["url"], CONFIG_ELN["token"])

# def get_current_time() -> str:
#     """
#     Get the current time formatted nicely for display.
#     """
#     # Get current datetime
#     now = datetime.now(ZoneInfo("Europe/Zurich"))
#     # Format as a human-readable sentence
#     formatted = now.strftime("%A, %B %d, %Y at %H:%M")
#     return f"It is currently {formatted}."

# @lru_cache(maxsize=5000)
# def _get_openbis_object(permId):
#     try:
#         obj = OPENBIS_SESSION.get_object(permId)
#         return obj
#     except ValueError as e:
#         return None

# @lru_cache(maxsize=5000)
# def _get_openbis_property_type(key):
#     try:
#         prop = OPENBIS_SESSION.get_property_type(key)
#         return prop
#     except ValueError as e:
#         return None

# @lru_cache(maxsize=100)
# def _get_openbis_objects_by_type(type):
#     openbis_objs = OPENBIS_SESSION.get_objects(type = type)
#     objs = {}
#     for openbis_obj in openbis_objs:
#         objs[openbis_obj.permId] = openbis_obj
#     return objs

# @lru_cache(maxsize=1000)
# def _get_openbis_experiments():
#     experiments = {}
#     openbis_experiments = OPENBIS_SESSION.get_experiments(type = "EXPERIMENT")
#     for exp in openbis_experiments:
#         experiments[exp.permId] = exp
#     return experiments

# @lru_cache(maxsize=1000)
# def _get_openbis_experiment(permId):
#     return OPENBIS_SESSION.get_experiment(permId)

# @lru_cache(maxsize=1000)
# def _get_substance_data(permId):
#     substance_data = {}
#     obj = _get_openbis_object(permId)
#     if obj:
#         substance_properties = obj.props.all()
#         for key, value in substance_properties.items():
#             if value:
#                 prop_type = _get_openbis_property_type(key)
#                 prop_datatype = prop_type.dataType
#                 if prop_datatype == "SAMPLE":
#                     if key == "molecules":
#                         molecules_data = []
#                         for molecule in value:
#                             molecule = _get_openbis_object(molecule)
#                             if molecule:
#                                 molecule_data = {}
#                                 molecule_properties = molecule.props.all()
#                                 for molecule_key, molecule_value in molecule_properties.items():
#                                     if molecule_value:
#                                         molecule_data[molecule_key] = molecule_value
                                
#                                 molecules_data.append(molecule_data)
                        
#                         substance_data["molecules"] = molecules_data
#                     else:
#                         if isinstance(value, list):
#                             prop_data = []
#                             for obj_permId in value:
#                                 openbis_object = _get_openbis_object(obj_permId)
#                                 if openbis_object:
#                                     openbis_object_properties = openbis_object.props.all()
#                                     openbis_name = openbis_object_properties.get("name", "Unknown")
#                                     prop_data.append({"name": openbis_name})
                            
#                             substance_data[key] = prop_data
                
#                 else:
#                     substance_data[key] = value
    
#     return substance_data

# @lru_cache(maxsize=100)
# def _get_experiment_data(permId):
#     experiment_data = {
#         "preparations": [],
#         "measurements": []
#     }
#     experiment = _get_openbis_experiment(permId)
#     exp_objects = experiment.get_objects()
    
#     for obj in exp_objects:
#         obj = _get_openbis_object(obj.permId)
#         if obj.type == "PREPARATION":
#             preparation_data = {
#                 "name": "",
#                 "process_steps": []
#             }
#             preparation_properties = obj.props.all()
#             preparation_data["name"] = preparation_properties["name"]
#             preparation_children = obj.children
            
#             for child in preparation_children:
#                 child_obj = _get_openbis_object(child)
#                 if child_obj.type == "PROCESS_STEP":
#                     process_step_properties = child_obj.props.all()
#                     process_step_data = {
#                         "name": "",
#                         "actions": [],
#                         "observables": [],
#                         "instrument": "",
#                         "sample_in": "",
#                         "sample_out": ""
#                     }
#                     process_step_data["name"] = process_step_properties["name"]
                    
#                     # Get sample in name
#                     process_step_parents = child_obj.parents
#                     for parent in process_step_parents:
#                         parent_obj = _get_openbis_object(parent)
#                         if parent_obj.type == "SAMPLE":
#                             process_step_data["sample_in"] = parent_obj.props["name"]
#                             break
                    
#                     # Get sample out name
#                     process_step_children = child_obj.children
#                     for process_step_child in process_step_children:
#                         process_step_child_obj = _get_openbis_object(process_step_child)
#                         if process_step_child_obj.type == "SAMPLE":
#                             process_step_data["sample_out"] = process_step_child_obj.props["name"]
#                             break
                    
#                     for step_key, step_value in process_step_properties.items():
#                         if step_value:
#                             if step_key == "actions":
#                                 for action_permId in step_value:
#                                     action_data = {
#                                         "name": "",
#                                         "type": "",
#                                         "description": "",
#                                         "duration": "",
#                                         "target_temperature": "",
#                                         "substrate_temperature": "",
#                                         "substance": "",
#                                         "component": "",
#                                         "component_settings": "",
#                                         "dosing_gas": "",
#                                         "pressure": "",
#                                         "sputter_ion": "",
#                                         "current": "",
#                                         "angle": "",
#                                         "comments": ""
#                                     }
                                    
#                                     action_obj = _get_openbis_object(action_permId)
#                                     action_properties = action_obj.props.all()
#                                     action_data["name"] = action_properties.pop("name")
#                                     action_data["type"] = action_obj.type.code
                                    
#                                     for action_key, action_value in action_properties.items():
#                                         if action_value:
#                                             prop_type = _get_openbis_property_type(action_key)
#                                             prop_datatype = prop_type.dataType
#                                             if prop_datatype == "SAMPLE":
#                                                 prop_sample = _get_openbis_object(action_value)
#                                                 prop_sample_name = prop_sample.props["name"]
#                                                 action_data[action_key] = prop_sample_name
#                                             else:
#                                                 action_data[action_key] = action_value
                                                
#                                     process_step_data["actions"].append(action_data)
                                    
#                             elif step_key == "observables":
#                                 for observable_permId in step_value:
#                                     observable_data = {
#                                         "name": "",
#                                         "type": "",
#                                         "description": "",
#                                         "duration": "",
#                                         "channel_name": "",
#                                         "component": "",
#                                         "component_settings": "",
#                                         "comments": ""
#                                     }
                                    
#                                     observable_obj = _get_openbis_object(observable_permId)
#                                     observable_properties = observable_obj.props.all()
#                                     observable_data["name"] = observable_properties.pop("name")
#                                     observable_data["type"] = observable_obj.type.code
                                    
#                                     for observable_key, observable_value in observable_properties.items():
#                                         if observable_value:
#                                             prop_type = _get_openbis_property_type(observable_key)
#                                             prop_datatype = prop_type.dataType
#                                             if prop_datatype == "SAMPLE":
#                                                 prop_sample = _get_openbis_object(observable_value)
#                                                 prop_sample_name = prop_sample.props["name"]
#                                                 observable_data[observable_key] = prop_sample_name
#                                             else:
#                                                 observable_data[observable_key] = observable_value
                                                
#                                     process_step_data["observables"].append(observable_data)
#                             elif step_key == "instrument":
#                                 instrument_obj = _get_openbis_object(step_value)
#                                 process_step_data["instrument"] = instrument_obj.props["name"]
#                             else:
#                                 process_step_data[step_key] = step_value

#                     preparation_data["process_steps"].append(process_step_data)
            
#             experiment_data["preparations"].append(preparation_data)
            
#         elif obj.type == "MEASUREMENT_SESSION":
#             pass
    
#     return experiment_data

# @lru_cache(maxsize=100)
# def _get_simulation_data(permId):
#     simulations_data = {
#         "simulations": [],
#         "measurements": []
#     }
#     experiment = _get_openbis_experiment(permId)
#     exp_objects = experiment.get_objects()
    
#     for obj in exp_objects:
#         obj = _get_openbis_object(obj.permId)
#         if obj.type in ["BAND_STRUCTURE", "GEOMETRY_OPTIMISATION", "VIBRATIONAL_SPECTROSCOPY", "UNCLASSIFIED_SIMULATION", "PDOS"]:
#             simulation_data = {
#                 "name": "",
#                 "wfms_uuid": "",
#                 "band_gap": "",
#                 "level_theory_method": "",
#                 "level_theory_parameters": "",
#                 "input_parameters": "",
#                 "output_parameters": "",
#                 "cell_optimisation_contraints": "",
#                 "cell_optimised": "",
#                 "driver_code": "",
#                 "constrained": "",
#                 "force_convergence_threshold": "",
#                 "description": "",
#                 "comments": "",
#                 "codes": [],
#                 "atomistic_model": []
#             }
            
#             simulation_properties = obj.props.all()
#             simulation_data["name"] = simulation_properties.pop("name")
#             simulation_data["type"] = obj.type.code
            
#             for sim_prop_key, sim_prop_value in simulation_properties.items():
#                 if sim_prop_value:
#                     sim_prop_type = _get_openbis_property_type(sim_prop_key)
#                     sim_prop_datatype = sim_prop_type.dataType
#                     if sim_prop_datatype == "SAMPLE":
#                         if isinstance(sim_prop_value, list):
#                             sim_prop_obj_data = []
#                             for sim_prop_obj in sim_prop_value:
#                                 sim_prop_obj = _get_openbis_object(sim_prop_obj)
#                                 if sim_prop_obj:
#                                     sim_prop_obj_properties = sim_prop_obj.props.all()
#                                     openbis_name = sim_prop_obj_properties.get("name", "Unknown")
#                                     sim_prop_obj_data.append({"name": openbis_name})
                            
#                             simulation_data[sim_prop_key] = sim_prop_obj_data
#                         else:
#                             prop_sample = _get_openbis_object(sim_prop_value)
#                             prop_sample_name = prop_sample.props["name"]
#                             simulation_data[sim_prop_key] = prop_sample_name
#                     else:
#                         simulation_data[sim_prop_key] = sim_prop_value

#             simulation_parents = obj.parents
#             for parent in simulation_parents:
#                 parent_obj = _get_openbis_object(parent)
#                 if parent_obj.type == "ATOMISTIC_MODEL":
#                     atomistic_model_properties = parent_obj.props.all()
#                     atomistic_model_data = {
#                         "molecules": [],
#                         "crystal_concepts": []
#                     }
                    
#                     for atom_model_key, atom_model_value in atomistic_model_properties.items():
#                         atomistic_model_data[atom_model_key] = atom_model_value
                    
#                     atomistic_model_parents = parent_obj.parents
#                     for atom_model_parent in atomistic_model_parents:
#                         atom_model_parent = _get_openbis_object(atom_model_parent)
                        
#                         if atom_model_parent.type == "MOLECULE":
#                             molecule_name = atom_model_parent.props["name"]
#                             atomistic_model_data["molecules"].append(molecule_name)
                            
#                         elif atom_model_parent.type == "CRYSTAL_CONCEPT":
#                             crystal_concept_name = atom_model_parent.props["name"]
#                             atomistic_model_data["crystal_concepts"].append(crystal_concept_name)
                            
#                     simulation_data["atomistic_model"].append(atomistic_model_data)
#                     break
                        
#             simulations_data["simulations"].append(simulation_data)
            
#         elif obj.type == "MEASUREMENT_SESSION":
#             pass
    
#     return simulations_data

# # Converts a Python function into a LangChain-compatible tool 
# # that the agent can call automatically
# @tool("get_substances")
# def _get_substances() -> List[Dict]:
#     """Returns information about all substances or molecules available in the inventory of openBIS."""
#     substances_data = []
#     substances_objects = _get_openbis_objects_by_type("SUBSTANCE")
#     for obj_permId, obj in substances_objects.items():
#         substance_data = _get_substance_data(obj_permId)
#         substances_data.append(substance_data)

#     return substances_data

# @tool("get_substances_by_name")
# def _get_substances_by_name(name: str) -> List[Dict]:
#     """Returns substance or molecule metadata available in the inventory of openBIS by name."""
#     substances_objects = _get_openbis_objects_by_type("SUBSTANCE")
#     substances_data = []
#     for obj_permId, obj in substances_objects.items():
#         substance_properties = obj.props.all()
#         if substance_properties["name"] == name:
#             substance_data = _get_substance_data(obj_permId)
#             substances_data.append(substance_data)

#     return substances_data
     
# @tool("get_substances_by_smiles")
# def _get_substances_by_smiles(smiles: str) -> List[Dict]:
#     """Returns substance or molecule metadata available in the inventory of openBIS by SMILES string."""
#     substances_data = []
#     substances_objects = _get_openbis_objects_by_type("SUBSTANCE")
#     for substance_obj_permId, substance_obj in substances_objects.items():
#         smiles_found = False
#         substance_properties = substance_obj.props.all()
#         for substance_prop_key, substance_prop_value in substance_properties.items():
#             if substance_prop_key == "molecules" and substance_prop_value:
#                 for molecule_permId in substance_prop_value:
#                     molecule_obj = _get_openbis_object(molecule_permId)
#                     if molecule_obj:
#                         molecule_properties = molecule_obj.props.all()
#                         if molecule_properties["smiles"] == smiles:
#                             smiles_found = True
                
#                 if smiles_found:
#                     substance_data = _get_substance_data(substance_obj_permId)
#                     substances_data.append(substance_data)

#     return substances_data

# @tool("get_substances_by_empa_identifier")
# def _get_substances_by_empa_identifier(empa_number: str, batch: str) -> List[Dict]:
#     """Returns substance or molecule metadata available in the inventory of openBIS by Empa number and batch, e.g., molecule 125a or 700b."""
#     substances_objects = _get_openbis_objects_by_type("SUBSTANCE")
#     substances_data = []
#     for obj_permId, obj in substances_objects.items():
#         substance_properties = obj.props.all()
#         if substance_properties["empa_number"] == empa_number and substance_properties["batch"] == batch:
#             substance_data = _get_substance_data(obj_permId)
#             substances_data.append(substance_data)
#             break

#     return substances_data

@tool("get_openbis_schema")
def _get_openbis_schema() -> dict:
    """
    Returns a dictionary that maps each OpenBIS object type to its corresponding list of property codes (attributes).
    This information is essential for understanding which filters can be applied to each object type and must always be 
    retrieved before querying OpenBIS objects to ensure valid filtering and accurate interpretation.
    """
    schema = {}
    try:
        # Get all object types
        object_types = {obj_type.code: obj_type for obj_type in OPENBIS_SESSION.get_object_types()}
        for type_code, obj_type in object_types.items():
            try:
                property_df = obj_type.get_property_assignments().df
                property_codes = [code.lower() for code in property_df["code"].tolist()]
                schema[type_code] = property_codes
            except Exception as e:
                schema[type_code] = f"Error fetching properties: {str(e)}"
    except Exception as e:
        return {"error": f"Failed to retrieve OpenBIS schema: {str(e)}"}

    return schema

@tool("get_openbis_objects")
def _get_openbis_objects(filters: dict = None, permId: str = None, type: str = None) -> list:
    """Return openBIS object data following this logic:
    If a permId is provided, retrieve the object directly by permId and ignore all filters and type.
    A permId typically has the format YYYYMMDDHHMMSSmmm-#### (e.g., 20250717103412473-3103).
    If a property value matches the permId format, try to resolve and return the linked object instead of the raw ID.
    If no permId is provided:
    Use the specified object type, which must be valid according to the openBIS schema.
    Apply the given filters based on the properties defined for that type.
    Only use properties that are assigned to the selected object type.
    If a filter references a property in permId format, resolve and return the referenced object, not just the ID."""
    
    if permId:
        openbis_object = OPENBIS_SESSION.get_object(permId)
        obj_properties = openbis_object.props.all()
        obj_properties["type"] = openbis_object.type.code
        obj_parents = openbis_object.get_parents()
        if obj_parents:
            obj_parents = list(obj_parents.df.permId)
        obj_properties["parents"] = obj_parents
        
        obj_children = openbis_object.get_children()
        if obj_children:
            obj_children = list(obj_children.df.permId)
        obj_properties["children"] = obj_children    
        list_objects = [obj_properties]
    else:
        if type:
            openbis_objects = OPENBIS_SESSION.get_objects(type = type)
        else:
            openbis_objects = OPENBIS_SESSION.get_objects()
        
        list_objects = []
        for obj in openbis_objects:
            obj_properties = obj.props.all()
            obj_properties["type"] = obj.type.code
            
            obj_parents = obj.get_parents()
            if obj_parents:
                obj_parents = list(obj_parents.df.permId)
            obj_properties["parents"] = obj_parents
            
            obj_children = obj.get_children()
            if obj_children:
                obj_children = list(obj_children.df.permId)
            obj_properties["children"] = obj_children
            
            if filters:
                match = all(k in obj_properties and obj_properties[k] == v for k, v in filters.items())
                if match and obj_properties not in list_objects:
                    list_objects.append(obj_properties)
            else:
                list_objects.append(obj_properties)
    
    return list_objects

@tool("get_experiments")
def _get_experiments() -> List[Dict]:
    """Returns a list of experiments, each with all the preparations, process steps, actions, observables and materials used.
    The materials include crystals, 2D materials, substances (molecules), among others. The types are returned in capital letters and
    with underscore instead of space, e.g. LIGHT_IRRADIATION. When ingesting the data convert them to normal words, e.g. from LIGHT_IRRADIATION to Light Irradiation."""
    
    experiments = _get_openbis_experiments()
    experiments_data = []
    for exp_permId, exp in experiments.items():
        experiment_data = _get_experiment_data(exp_permId)
        experiments_data.append(experiment_data)
    
    return experiments_data

@tool("get_simulations")
def _get_simulations() -> List[Dict]:
    """Returns a list of simulations, each with all the atomistic models, aiida nodes, codes, used materials (molecules, crystal concepts, 2D materials, among others).
    There are different simulations types such as band structure, geometry optimisation, vibrational spectroscopy, among others. The types are returned in capital letters and
    with underscore instead of space, e.g. BAND_STRUCTURE. When ingesting the data convert them to normal words, e.g. from BAND_STRUCTURE to Band Structure."""
    
    experiments = _get_openbis_experiments()
    simulations_data = []
    for exp_permId, exp in experiments.items():
        simulation_data = _get_simulation_data(exp_permId)
        simulations_data.append(simulation_data)
    
    return simulations_data

from langchain.agents import tool

class OpenBISAgent():
    """
    An agent that can answer questions about openBIS inventories.
    It uses a LangChain model and tools to interact with the openBIS session.
    """
    def __init__(self, google_api_key):
        self.google_api_key = google_api_key
        self.llm_model = ChatGoogleGenerativeAI(model = "models/gemini-2.5-flash-lite-preview-06-17", google_api_key = self.google_api_key)
        self.system_prompt = """
            You are an helpful assistant that can answer questions about experiments, simulations, and all the inventories inside 
            openBIS. In the case you are asked for an object that contains links to other objects, get all the information
            about those objects instead of giving the identifier. {get_current_time()}
        """
        
        self._tools = [
            # _get_substances,
            # _get_substances_by_name,
            # _get_substances_by_smiles,
            # _get_substances_by_empa_identifier,
            # _get_crystals,
            # _get_experiments,
            # _get_simulations,
            _get_openbis_schema,
            _get_openbis_objects
        ]
        
        self.ai_agent = create_react_agent(
            model = self.llm_model,
            tools = self._tools,
            checkpointer = InMemorySaver(),
            prompt = self.system_prompt
        )
        
        self.agent_config = {"configurable": {"thread_id": "1"}}

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
        response = self.ai_agent.invoke(
            {
                "messages": [
                    {"role": "user", "content": user_prompt}
                ]
            },
            self.agent_config
        )

        return response["messages"]