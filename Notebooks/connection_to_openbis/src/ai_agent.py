import src.utils as utils
import asyncio
from zoneinfo import ZoneInfo
from datetime import datetime
from functools import lru_cache
from langchain_core.tools import tool
from langgraph.checkpoint.memory import InMemorySaver
from langgraph.prebuilt import ToolNode, tools_condition, create_react_agent
from langchain_google_genai import ChatGoogleGenerativeAI
from typing import List, Dict, Literal, Annotated
from typing_extensions import TypedDict
from langgraph.graph import StateGraph, START, END
from langgraph.graph.message import add_messages
from uuid import uuid4

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

@lru_cache(maxsize=100)
def _get_experiment_data(permId):
    experiment_data = {
        "preparations": [],
        "measurements": []
    }
    experiment = _get_openbis_experiment(permId)
    exp_objects = experiment.get_objects()
    
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
    
    return experiment_data

@lru_cache(maxsize=100)
def _get_simulation_data(permId):
    simulations_data = {
        "simulations": [],
        "measurements": []
    }
    experiment = _get_openbis_experiment(permId)
    exp_objects = experiment.get_objects()
    
    for obj in exp_objects:
        obj = _get_openbis_object(obj.permId)
        if obj.type in ["BAND_STRUCTURE", "GEOMETRY_OPTIMISATION", "VIBRATIONAL_SPECTROSCOPY", "UNCLASSIFIED_SIMULATION", "PDOS"]:
            simulation_data = {
                "name": "",
                "wfms_uuid": "",
                "band_gap": "",
                "level_theory_method": "",
                "level_theory_parameters": "",
                "input_parameters": "",
                "output_parameters": "",
                "cell_optimisation_contraints": "",
                "cell_optimised": "",
                "driver_code": "",
                "constrained": "",
                "force_convergence_threshold": "",
                "description": "",
                "comments": "",
                "codes": [],
                "atomistic_model": []
            }
            
            simulation_properties = obj.props.all()
            simulation_data["name"] = simulation_properties.pop("name")
            simulation_data["type"] = obj.type.code
            
            for sim_prop_key, sim_prop_value in simulation_properties.items():
                if sim_prop_value:
                    sim_prop_type = _get_openbis_property_type(sim_prop_key)
                    sim_prop_datatype = sim_prop_type.dataType
                    if sim_prop_datatype == "SAMPLE":
                        if isinstance(sim_prop_value, list):
                            sim_prop_obj_data = []
                            for sim_prop_obj in sim_prop_value:
                                sim_prop_obj = _get_openbis_object(sim_prop_obj)
                                if sim_prop_obj:
                                    sim_prop_obj_properties = sim_prop_obj.props.all()
                                    openbis_name = sim_prop_obj_properties.get("name", "Unknown")
                                    sim_prop_obj_data.append({"name": openbis_name})
                            
                            simulation_data[sim_prop_key] = sim_prop_obj_data
                        else:
                            prop_sample = _get_openbis_object(sim_prop_value)
                            prop_sample_name = prop_sample.props["name"]
                            simulation_data[sim_prop_key] = prop_sample_name
                    else:
                        simulation_data[sim_prop_key] = sim_prop_value

            simulation_parents = obj.parents
            for parent in simulation_parents:
                parent_obj = _get_openbis_object(parent)
                if parent_obj.type == "ATOMISTIC_MODEL":
                    atomistic_model_properties = parent_obj.props.all()
                    atomistic_model_data = {
                        "molecules": [],
                        "crystal_concepts": []
                    }
                    
                    for atom_model_key, atom_model_value in atomistic_model_properties.items():
                        atomistic_model_data[atom_model_key] = atom_model_value
                    
                    atomistic_model_parents = parent_obj.parents
                    for atom_model_parent in atomistic_model_parents:
                        atom_model_parent = _get_openbis_object(atom_model_parent)
                        
                        if atom_model_parent.type == "MOLECULE":
                            molecule_name = atom_model_parent.props["name"]
                            atomistic_model_data["molecules"].append(molecule_name)
                            
                        elif atom_model_parent.type == "CRYSTAL_CONCEPT":
                            crystal_concept_name = atom_model_parent.props["name"]
                            atomistic_model_data["crystal_concepts"].append(crystal_concept_name)
                            
                    simulation_data["atomistic_model"].append(atomistic_model_data)
                    break
                        
            simulations_data["simulations"].append(simulation_data)
            
        elif obj.type == "MEASUREMENT_SESSION":
            pass
    
    return simulations_data

# Converts a Python function into a LangChain-compatible tool 
# that the agent can call automatically
@tool("query_openbis_objects")
def _query_openbis_objects(permId: str = None, type: str = None, filters: dict = None) -> List[Dict]:
    """
    Returns list of objects obtained from openBIS. 
    
    This tool can return an object by permanent identifier (permID) or by type.
    
    When permID is given, it searches for the object and returns its properties. When it is not given but instead the type
    is given, it searches for all the objects of that type and returns a list with them. In the case of returning objects by type, 
    if filters are given, objects should be filtered in order to only return the ones that have certain attributes.
    
    Examples: 
        - User: Find me all molecules inside openBIS -> permId = None and type = MOLECULE
        - User: Give me a small summary of object 20250717114928650-3593 -> permId = 20250717114928650-3593 and type = None
        - User: Get me molecule 704 -> type = MOLECULE and filters = {'name': 704}
        - User: Get me substance 704a -> type = SUBSTANCE and filters = {'empa_number': 704, 'batch': 'a'}
        - User: Get me all samples inside openBIS -> type = SAMPLE
    
    - type should be selected among the types found in the get_openbis_schema tool
    
    - filters should be considered according to the type. It means that if a certain type is selected, the filters can only 
    be properties that belong to that type. Example: MOLECULE contains name, smiles, sum_formula, cas_number, iupac_name, and comments.
    It means that description cannot be used. 
        
    - permId allways follows this structure: YYYYMMDDHHMMSSmmm-####
    """
    
    if permId:
        openbis_object = _get_openbis_object(permId)
        if openbis_object:
            openbis_object_data = openbis_object.props.all()
            return [openbis_object_data]
        else:
            return f"Cannot find the object with permId {permId}."
    else:
        objects = _get_openbis_objects_by_type(type)
        list_objects = []
        for obj_permId, obj in objects.items():
            obj_props = obj.props.all()
            if filters:
                match = all(k in obj_props and obj_props[k] == v for k, v in filters.items())
                if match:
                    list_objects.append(obj_props)
            else:
                list_objects.append(obj_props)
        
        if list_objects:
            return list_objects
        else:
            return f"Cannot find objects with type {type} and filters {filters}."

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

@tool("get_openbis_schema")
def _get_openbis_schema() -> dict:
    """
    Returns a dictionary where keys are OpenBIS object types and values are lists of property codes (filters/attributes).
    This schema should always be fetched before filtering, so the model knows which filters are valid for each object type.
    
    Example output:
    {
        "MOLECULE": ["name", "smiles", "sum_formula"],
        "SUBSTANCE": ["empa_number", "batch"],
        ...
    }
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

# @tool("get_openbis_objects")
# def _get_openbis_objects(filters: dict = None, permId: str = None, type: str = None) -> list:
#     """Return openBIS object data following this logic:
#     If a permId is provided, retrieve the object directly by permId and ignore all filters and type.
#     A permId typically has the format YYYYMMDDHHMMSSmmm-#### (e.g., 20250717103412473-3103).
#     If a property value matches the permId format, try to resolve and return the linked object instead of the raw ID.
#     If no permId is provided:
#     Use the specified object type, which must be valid according to the openBIS schema.
#     Apply the given filters based on the properties defined for that type.
#     Only use properties that are assigned to the selected object type.
#     If a filter references a property in permId format, resolve and return the referenced object, not just the ID."""
    
#     if permId:
#         openbis_object = OPENBIS_SESSION.get_object(permId)
#         obj_properties = openbis_object.props.all()
#         obj_properties["type"] = openbis_object.type.code
#         obj_parents = openbis_object.get_parents()
#         if obj_parents:
#             obj_parents = list(obj_parents.df.permId)
#         obj_properties["parents"] = obj_parents
        
#         obj_children = openbis_object.get_children()
#         if obj_children:
#             obj_children = list(obj_children.df.permId)
#         obj_properties["children"] = obj_children    
#         list_objects = [obj_properties]
#     else:
#         if type:
#             openbis_objects = OPENBIS_SESSION.get_objects(type = type)
#         else:
#             openbis_objects = OPENBIS_SESSION.get_objects()
        
#         list_objects = []
#         for obj in openbis_objects:
#             obj_properties = obj.props.all()
#             obj_properties["type"] = obj.type.code
            
#             obj_parents = obj.get_parents()
#             if obj_parents:
#                 obj_parents = list(obj_parents.df.permId)
#             obj_properties["parents"] = obj_parents
            
#             obj_children = obj.get_children()
#             if obj_children:
#                 obj_children = list(obj_children.df.permId)
#             obj_properties["children"] = obj_children
            
#             if filters:
#                 match = all(k in obj_properties and obj_properties[k] == v for k, v in filters.items())
#                 if match and obj_properties not in list_objects:
#                     list_objects.append(obj_properties)
#             else:
#                 list_objects.append(obj_properties)
    
#     return list_objects

class State(TypedDict):
    messages: Annotated[list, add_messages]

@tool("get_openbis_schema")
def _get_openbis_schema() -> dict:
    """
    Returns a dictionary where keys are OpenBIS object types and values are lists of property codes (filters/attributes).
    This schema should always be fetched before calling the tool for accessing object data and filter object data (query_openbis_objects tool), 
    so the model knows which filters are valid for each object type.
    
    Example output:
    {
        "MOLECULE": ["name", "smiles", "sum_formula"],
        "SUBSTANCE": ["empa_number", "batch"],
        ...
    }
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

class OpenBISAgent():
    """
    An agent that can answer questions about openBIS inventories.
    It uses a LangChain model and tools to interact with the openBIS session.
    """
    def __init__(self, google_api_key):
        self.google_api_key = google_api_key
        self.llm_model = ChatGoogleGenerativeAI(model = "models/gemini-2.5-flash", google_api_key = self.google_api_key)
        
        self.system_prompt = f"""
            You are a helpful assistant that can answer questions about experiments, simulations, 
            and all the inventories inside openBIS. 
            In the case you are asked for an object that contains links to other objects, get all the information
            about those objects instead of giving the identifier. 
            Today is {get_current_time()}.
        """
        
        self._tools = [
            _get_experiments,
            _get_simulations,
            _query_openbis_objects,
            _get_openbis_schema
        ]
        
        llm_with_tools = self.llm_model.bind_tools(self._tools)
        def chatbot(state: State):
            return {"messages": [llm_with_tools.invoke(state["messages"])]}
        
        # Define the state machine
        # Add the nodes
        graph_builder = StateGraph(State)
        graph_builder.add_node("chatbot", chatbot)
        tool_node = ToolNode(tools = self._tools)
        graph_builder.add_node("tools", tool_node)
        graph_builder.add_conditional_edges("chatbot", tools_condition)
        
        # Build the edges. Make sure that the bot always gets the openBIS schema before using the other tools
        graph_builder.add_edge(START, "chatbot")
        graph_builder.add_edge("tools", "chatbot")
        
        memory = InMemorySaver()
        self.graph = graph_builder.compile(checkpointer=memory)
        self.agent_config = {"configurable": {"thread_id": uuid4()}}

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
        events = self.graph.stream(
            {"messages": [{"role": "user", "content": user_prompt}]},
            self.agent_config,
            stream_mode="values",
        )
        response = []
        for event in events:
            response.append(event["messages"][-1])
        return response

    def get_graph_replay(self):
        to_replay = None
        for state in self.graph.get_state_history(self.agent_config):
            print("Num Messages: ", len(state.values["messages"]), "Next: ", state.next)
            print("-" * 80)
            if len(state.values["messages"]) == 10:
                # We are somewhat arbitrarily selecting a specific state based on the number of chat messages in the state.
                to_replay = state
        
        return to_replay

# class OpenBISAgent():
#     """
#     An agent that can answer questions about openBIS inventories.
#     It uses a LangChain model and tools to interact with the openBIS session.
#     """
#     def __init__(self, google_api_key):
#         self.google_api_key = google_api_key
#         self.llm_model = ChatGoogleGenerativeAI(model = "models/gemini-2.5-flash", google_api_key = self.google_api_key)
#         self.system_prompt = """
#             You are an helpful assistant that can answer questions about experiments, simulations, and all the inventories inside 
#             openBIS. In the case you are asked for an object that contains links to other objects, get all the information
#             about those objects instead of giving the identifier. {get_current_time()}
#         """
        
#         self._tools = [
#             _get_experiments,
#             _get_simulations,
#             _query_openbis_objects,
#             _get_openbis_schema
#         ]
        
#         self.ai_agent = create_react_agent(
#             model = self.llm_model,
#             tools = self._tools,
#             checkpointer = InMemorySaver(),
#             prompt = self.system_prompt
#         )
        
#         self.agent_config = {"configurable": {"thread_id": uuid4()}}

#     def ask_question(self, user_prompt: str):
#         """
#         Ask a question to the agent and get a response.
        
#         Args:
#             user_prompt (str): The question to ask the agent.
        
#         Returns:
#             str: The response from the agent.
#         """
#         # Builds a LangGraph agent that interleaves reasoning and tool execution (ReAct pattern)
#         # Runs the full workflow, returning a state object 
#         # where the last message contains the model’s response after any tool invocations
#         response = self.ai_agent.invoke(
#             {
#                 "messages": [
#                     {"role": "user", "content": user_prompt}
#                 ]
#             },
#             self.agent_config
#         )

#         return response["messages"]