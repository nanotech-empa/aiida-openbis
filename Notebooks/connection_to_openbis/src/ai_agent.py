import src.utils as utils
import asyncio
from zoneinfo import ZoneInfo
from datetime import datetime
from functools import lru_cache
from langchain_core.tools import tool
from langgraph.checkpoint.memory import InMemorySaver
from langgraph.prebuilt import create_react_agent
from langchain_google_genai import ChatGoogleGenerativeAI

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

# Converts a Python function into a LangChain-compatible tool 
# that the agent can call automatically
@tool("substances")
def _get_substances():
    """Call to get information about substances or molecules available in the inventory of openBIS. 
    This function gets all the substances or molecules available in openBIS and its useful to answer questions about all of them or just a specific one."""
    substances_objects = _get_openbis_objects_by_type("SUBSTANCE")
    substances_info = ""
    for obj_permId, obj in substances_objects.items():
        substance_properties = obj.props.all()
        for key, value in substance_properties.items():
            if value:
                prop_type = _get_openbis_property_type(key)
                prop_label = prop_type.label
                prop_datatype = prop_type.dataType
                if prop_datatype == "SAMPLE":
                    if key == "molecules":
                        molecule_info = "Molecules: "
                        for molecule in value:
                            molecule = _get_openbis_object(molecule)
                            if molecule:
                                molecule_properties = molecule.props.all()
                                for molecule_key, molecule_value in molecule_properties.items():
                                    if molecule_value:
                                        molecule_prop_type = _get_openbis_property_type(key)
                                        molecule_prop_label = molecule_prop_type.label
                                        molecule_info += f"{molecule_prop_label}: {molecule_value}\n"
                        
                        value = molecule_info
                    else:
                        if isinstance(value, list):
                            for obj_permId in value:
                                openbis_object = _get_openbis_object(obj_permId)
                                if openbis_object:
                                    openbis_object_properties = openbis_object.props.all()
                                    openbis_name = openbis_object_properties.get("name", "Unknown")
                                    value += f"Name: {openbis_name}\n"

                substances_info += f"{prop_label}: {value}\n"

    return substances_info

@tool("crystals")
def _get_crystals():
    """Call to get information about crystals available in the inventory of openBIS. 
    This function gets all the crystals available in openBIS and its useful to answer questions about all of them or just a specific one."""
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

@tool("experiments")
def _get_experiments():
    """Call to get information about experiments stored in openBIS. Experiments concern everything that was
    performed in the labs, like sample preparation (annealing, sputtering, deposition, dosing, etc.)."""
    experiments = _get_openbis_experiments()
    all_experiments_info = ""
    for exp_permId, exp in experiments.items():
        exp_objects = exp.get_objects()
        experiment_info = ""
        for obj in exp_objects:
            obj = _get_openbis_object(obj.permId)
            if obj.type == "PREPARATION":
                preparation_properties = obj.props.all()
                preparation_info = "Preparation: "
                for key, value in preparation_properties.items():
                    if value:
                        prop_type = _get_openbis_property_type(key)
                        prop_label = prop_type.label
                        preparation_info += f"{prop_label}: {value}\n"
                preparation_children = obj.children
                for child in preparation_children:
                    child_obj = _get_openbis_object(child)
                    if child_obj.type == "PROCESS_STEP":
                        step_properties = child_obj.props.all()
                        step_info = "Process step: "
                        for key, value in step_properties.items():
                            if value:
                                if key == "actions":
                                    actions_objs = []
                                    actions_info = "Actions: "
                                    for action_permId in value:
                                        action_obj = _get_openbis_object(action_permId)
                                        action_type = action_obj.type
                                        action_properties = action_obj.props.all()
                                        actions_info += f"{action_type}: "
                                        for key, value in action_properties.items():
                                            if value:
                                                if key == "substance":
                                                    substance_info = "Substance: "
                                                    prop_type = _get_openbis_property_type(key)
                                                    prop_label = prop_type.label
                                                    prop_datatype = prop_type.dataType
                                                    if prop_datatype == "SAMPLE":
                                                        if key == "molecules":
                                                            molecule_info = "Molecules: "
                                                            for molecule in value:
                                                                molecule = _get_openbis_object(molecule)
                                                                if molecule:
                                                                    molecule_properties = molecule.props.all()
                                                                    for molecule_key, molecule_value in molecule_properties.items():
                                                                        if molecule_value:
                                                                            molecule_prop_type = _get_openbis_property_type(key)
                                                                            molecule_prop_label = molecule_prop_type.label
                                                                            molecule_info += f"{molecule_prop_label}: {molecule_value}\n"
                                                                            
                                                            substance_info += molecule_info
                                                
                                                elif key == "component":
                                                    component_obj = _get_openbis_object(value)
                                                    component_properties = component_obj.props.all()
                                                    component_name = component_properties.get("name", "Unknown")
                                                    step_info += f"Component: {component_name}\n"
                                                    
                                                else:
                                                    prop_type = _get_openbis_property_type(key)
                                                    prop_label = prop_type.label
                                                    actions_info += f"{prop_label}: {value}\n"
                                    step_info += actions_info
                                        
                                elif key == "observables":
                                    pass
                                elif key == "instrument":
                                    instrument_obj = _get_openbis_object(value)
                                    instrument_properties = instrument_obj.props.all()
                                    instrument_name = instrument_properties.get("name", "Unknown")
                                    step_info += f"Instrument: {instrument_name}\n"
                                else:
                                    prop_type = _get_openbis_property_type(key)
                                    prop_label = prop_type.label
                                    step_info += f"{prop_label}: {value}\n"
                
                    preparation_info += step_info
                
                experiment_info += preparation_info
                
            elif obj.type == "MEASUREMENT_SESSION":
                pass
        
        all_experiments_info += f"Experiment: {experiment_info}"
    
    return all_experiments_info

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