import src.utils as utils
from zoneinfo import ZoneInfo
from datetime import datetime
from functools import lru_cache
from langchain_core.tools import tool
from langgraph.checkpoint.memory import InMemorySaver
from langgraph.prebuilt import ToolNode, tools_condition
from langchain_google_genai import ChatGoogleGenerativeAI
from typing import List, Dict, Annotated
from typing_extensions import TypedDict
from langgraph.graph import StateGraph, START
from langgraph.graph.message import add_messages
from langchain_core.messages import SystemMessage
from uuid import uuid4

# OpenBIS configuration and session setup
CONFIG_ELN = utils.get_aiidalab_eln_config()
OPENBIS_SESSION, SESSION_DATA = utils.connect_openbis(
    CONFIG_ELN["url"], CONFIG_ELN["token"]
)


def get_current_time() -> str:
    """
    Get the current time formatted nicely for display.
    """
    # Get current datetime
    now = datetime.now(ZoneInfo("Europe/Zurich"))
    formatted = now.strftime("%A, %B %d, %Y at %H:%M")
    return f"It is currently {formatted}."


# Format as a human-readable sentence
@lru_cache(maxsize=5000)
def _get_openbis_object(permId):
    try:
        obj = OPENBIS_SESSION.get_object(permId)
        return obj
    except ValueError:
        return None


@lru_cache(maxsize=5000)
def _get_openbis_property_type(key):
    try:
        prop = OPENBIS_SESSION.get_property_type(key)
        return prop
    except ValueError:
        return None


@lru_cache(maxsize=100)
def _get_openbis_objects_by_type(type):
    openbis_objs = OPENBIS_SESSION.get_objects(type=type)
    objs = {}
    for openbis_obj in openbis_objs:
        objs[openbis_obj.permId] = openbis_obj
    return objs


@lru_cache(maxsize=1000)
def _get_openbis_experiments():
    experiments = {}
    openbis_experiments = OPENBIS_SESSION.get_experiments(type="EXPERIMENT")
    for exp in openbis_experiments:
        experiments[exp.permId] = exp
    return experiments


@lru_cache(maxsize=1000)
def _get_openbis_experiment(permId):
    return OPENBIS_SESSION.get_experiment(permId)


@lru_cache(maxsize=100)
def _get_experiment_data(permId):
    experiment_data = {"preparations": [], "measurements": []}
    experiment = _get_openbis_experiment(permId)
    exp_objects = experiment.get_objects()

    for obj in exp_objects:
        obj = _get_openbis_object(obj.permId)
        if obj.type == "PREPARATION":
            preparation_data = {"name": "", "process_steps": []}
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
                        "sample_out": "",
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
                            process_step_data["sample_out"] = (
                                process_step_child_obj.props["name"]
                            )
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
                                        "comments": "",
                                    }

                                    action_obj = _get_openbis_object(action_permId)
                                    action_properties = action_obj.props.all()
                                    action_data["name"] = action_properties.pop("name")
                                    action_data["type"] = action_obj.type.code

                                    for (
                                        action_key,
                                        action_value,
                                    ) in action_properties.items():
                                        if action_value:
                                            prop_type = _get_openbis_property_type(
                                                action_key
                                            )
                                            prop_datatype = prop_type.dataType
                                            if prop_datatype == "SAMPLE":
                                                prop_sample = _get_openbis_object(
                                                    action_value
                                                )
                                                prop_sample_name = prop_sample.props[
                                                    "name"
                                                ]
                                                action_data[action_key] = (
                                                    prop_sample_name
                                                )
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
                                        "comments": "",
                                    }

                                    observable_obj = _get_openbis_object(
                                        observable_permId
                                    )
                                    observable_properties = observable_obj.props.all()
                                    observable_data["name"] = observable_properties.pop(
                                        "name"
                                    )
                                    observable_data["type"] = observable_obj.type.code

                                    for (
                                        observable_key,
                                        observable_value,
                                    ) in observable_properties.items():
                                        if observable_value:
                                            prop_type = _get_openbis_property_type(
                                                observable_key
                                            )
                                            prop_datatype = prop_type.dataType
                                            if prop_datatype == "SAMPLE":
                                                prop_sample = _get_openbis_object(
                                                    observable_value
                                                )
                                                prop_sample_name = prop_sample.props[
                                                    "name"
                                                ]
                                                observable_data[observable_key] = (
                                                    prop_sample_name
                                                )
                                            else:
                                                observable_data[observable_key] = (
                                                    observable_value
                                                )

                                    process_step_data["observables"].append(
                                        observable_data
                                    )
                            elif step_key == "instrument":
                                instrument_obj = _get_openbis_object(step_value)
                                process_step_data["instrument"] = instrument_obj.props[
                                    "name"
                                ]
                            else:
                                process_step_data[step_key] = step_value

                    preparation_data["process_steps"].append(process_step_data)

            experiment_data["preparations"].append(preparation_data)

        elif obj.type == "MEASUREMENT_SESSION":
            pass

    return experiment_data


@lru_cache(maxsize=100)
def _get_simulation_data(permId):
    simulations_data = {"simulations": [], "measurements": []}
    experiment = _get_openbis_experiment(permId)
    exp_objects = experiment.get_objects()

    for obj in exp_objects:
        obj = _get_openbis_object(obj.permId)
        if obj.type in [
            "BAND_STRUCTURE",
            "GEOMETRY_OPTIMISATION",
            "VIBRATIONAL_SPECTROSCOPY",
            "UNCLASSIFIED_SIMULATION",
            "PDOS",
        ]:
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
                "atomistic_model": [],
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
                                    openbis_name = sim_prop_obj_properties.get(
                                        "name", "Unknown"
                                    )
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
                    atomistic_model_data = {"molecules": [], "crystal_concepts": []}

                    for (
                        atom_model_key,
                        atom_model_value,
                    ) in atomistic_model_properties.items():
                        atomistic_model_data[atom_model_key] = atom_model_value

                    atomistic_model_parents = parent_obj.parents
                    for atom_model_parent in atomistic_model_parents:
                        atom_model_parent = _get_openbis_object(atom_model_parent)

                        if atom_model_parent.type == "MOLECULE":
                            molecule_name = atom_model_parent.props["name"]
                            atomistic_model_data["molecules"].append(molecule_name)

                        elif atom_model_parent.type == "CRYSTAL_CONCEPT":
                            crystal_concept_name = atom_model_parent.props["name"]
                            atomistic_model_data["crystal_concepts"].append(
                                crystal_concept_name
                            )

                    simulation_data["atomistic_model"].append(atomistic_model_data)
                    break

            simulations_data["simulations"].append(simulation_data)

        elif obj.type == "MEASUREMENT_SESSION":
            pass

    return simulations_data


# Converts a Python function into a LangChain-compatible tool
# that the agent can call automatically
@tool("query_openbis_objects")
def query_openbis_objects(
    permId: str = None, type: str = None, filters: dict = None
) -> List[Dict]:
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
        - User: Get me substance with empa number 704 and batch a -> type = SUBSTANCE and filters = {'empa_number': 704, 'batch': 'a'}
        - User: Get me all samples inside openBIS -> type = SAMPLE

    Args:
        type (str): Should be selected among the types found in the openBIS objects schema. It is optional.

        filters (dict): Should be considered according to the type. It means that if a certain type is selected,
        the filters can only be properties that belong to that type. Example: MOLECULE contains name, smiles,
        sum_formula, cas_number, iupac_name, and comments. It means that description cannot be used. It is optional.

        permId (str): Always follows this pattern: YYYYMMDDHHMMSSmmm-####. It is optional

    Returns:
        list[dict]: A list of dictionaries that contains objects' data. Each dictionary contains the identifier, the
        type of the object, the registration date, the properties, the parent objects identifiers and the children objects
        identifiers.
        Example:
        [
            {
                "permId": "20250717115421455-9999"
                "type": "PROCESS_STEP",
                "registration_date": "2025-07-22 09:58:40",
                "properties": {
                    "name": "Deposition",
                    "actions": ["20250717115421436-23241", "20250717115421436-984"],
                    "observables": ["20250717115421436-2222", "20250717115421436-3333"],
                    "instrument": []
                },
                "parents_objects": ["20240517115421436-5432", "20230217115421436-3113"],
                "children_objects": ["20250517115421445-5432", "20250317115421438-231"]
            }
            ...
        ]

    """

    if permId:
        openbis_object = _get_openbis_object(permId)
        if openbis_object:
            obj_parents = openbis_object.get_parents()
            if obj_parents:
                obj_parents = list(obj_parents.df.permId)

            obj_children = openbis_object.get_children()
            if obj_children:
                obj_children = list(obj_children.df.permId)

            obj_props = openbis_object.props.all()

            obj_dict = {
                "type": type,
                "permId": permId,
                "registration_date": openbis_object.registrationDate,
                "properties": obj_props,
                "parent_objects": obj_parents,
                "children_objects": obj_children,
            }
            return [obj_dict]
        else:
            return f"Cannot find the object with permId {permId}."
    else:
        objects = _get_openbis_objects_by_type(type)
        list_objects = []
        for obj_permId, obj in objects.items():
            obj_props = obj.props.all()
            obj_dict = {
                "permId": obj.permId,
                "type": type,
                "registration_date": obj.registrationDate,
                "properties": obj_props,
            }

            if filters:
                match = True
                for k, v in filters.items():
                    if isinstance(v, int) or isinstance(v, float):
                        int_v = int(v)
                        if v == int_v:
                            v = str(int_v)
                        else:
                            v = str(v)

                    if k not in obj_props or obj_props[k] != v:
                        match = False
                        break

                if match:
                    list_objects.append(obj_dict)
            else:
                list_objects.append(obj_dict)

            if len(list_objects) == 10:
                return list_objects

        if list_objects:
            return list_objects
        else:
            return f"Cannot find objects with type {type} and filters {filters}."


@tool("query_openbis_experiments")
def query_openbis_experiments() -> List[Dict]:
    """
    Returns a list of experiments, each with the list of permanent identifiers of the objects
    that are saved within it. Inside the experiments, there could be sample preparations, simulations,
    measurements, substances, etc.

    Use this tool whenever you are asked about experiments, or when you need to find experiments
    related to a specific object (like a substance, sample, or measurement).

    Example:
        - Get me all the experiments that uses substance 704a.
        - Find me all the experiments in openBIS.
    """

    experiments = _get_openbis_experiments()
    experiments_data = []
    for exp_permId, exp in experiments.items():
        exp_name = exp.props["name"]
        exp_objects = exp.get_objects()
        objects_ids = []
        for obj in exp_objects:
            objects_ids.append(obj.permId)
        experiments_data.append({exp_name: objects_ids})

    return experiments_data


@tool("get_simulations")
def get_simulations() -> List[Dict]:
    """Returns a list of simulations, each with all the atomistic models, aiida nodes, codes, used materials (molecules, crystal concepts, 2D materials, among others).
    There are different simulations types such as band structure, geometry optimisation, vibrational spectroscopy, among others. The types are returned in capital letters and
    with underscore instead of space, e.g. BAND_STRUCTURE. When ingesting the data convert them to normal words, e.g. from BAND_STRUCTURE to Band Structure."""

    experiments = _get_openbis_experiments()
    simulations_data = []
    for exp_permId, exp in experiments.items():
        simulation_data = _get_simulation_data(exp_permId)
        simulations_data.append(simulation_data)

    return simulations_data


def get_openbis_schema() -> Dict:
    """
    Returns a dictionary that contains openBIS object types as keys and list of dictionaries as values. Every dictionary
    from this list contains two elements, the code of the property and the datatype of the property.
    This schema should be used, so the model knows which filters are valid for each object type and the datatypes of those
    filters.

    Example output:
    {
        "MOLECULE":
            [
                {
                    "code": "name",
                    "datatype": "varchar"
                },
                {
                    "code": "smiles",
                    "datatype": "varchar"
                },
            ],
        "SUBSTANCE":
            [
                {
                    "code": "name",
                    "datatype": "varchar"
                },
                {
                    "code": "empa_number",
                    "datatype": "integer"
                },
            ]
        ...
    }
    """
    schema = {}
    try:
        # Get all object types
        object_types = {
            obj_type.code: obj_type for obj_type in OPENBIS_SESSION.get_object_types()
        }
        for type_code, obj_type in object_types.items():
            try:
                properties_df = obj_type.get_property_assignments().df
                properties_codes = properties_df["code"]
                properties_datatypes = properties_df["dataType"]
                properties_info = []
                for code, datatype in zip(properties_codes, properties_datatypes):
                    property_info = {"code": code.lower(), "datatype": datatype.lower()}
                    properties_info.append(property_info)
                schema[type_code] = properties_info
            except Exception as e:
                schema[type_code] = f"Error fetching properties: {str(e)}"
    except Exception as e:
        return {"error": f"Failed to retrieve OpenBIS schema: {str(e)}"}

    return schema


class State(TypedDict):
    messages: Annotated[list, add_messages]


class OpenBISAgent:
    """
    An agent that can answer questions about openBIS inventories.
    It uses a LangChain model and tools to interact with the openBIS session.
    """

    def __init__(self, google_api_key):
        self.google_api_key = google_api_key
        self.llm_model = ChatGoogleGenerativeAI(
            model="models/gemini-2.5-flash", google_api_key=self.google_api_key
        )

        self.system_prompt = f"""
            You are a helpful assistant that can answer questions about experiments, simulations,
            and all the inventories inside openBIS.
            This is the schema containing the current openBIS object types and their properties/attributes:
            {get_openbis_schema()}. Use this schema whenever you are asked about objects in openBIS in order to understand
            their structure and types.

            In the case you are asked for an object that contains links to other objects, get all the information
            about those objects instead of giving the identifier. Example: You get a substance that contains a property named
            molecule which has the permID 20250717114928650-3593. Then, search for the object by permId and get its data.

            In case you are asked about an experiment that contains a certain object, first find the object details using
            query_openbis_objects tool and then query_openbis_experiments tool. This way you first get data about the object
            and then you get data about the experiment and you are able to combine both information to give an answer.
            Example:
                User: Is there substance with empa number 704 and batch a? If yes, find the experiments that use that and summarize them.
                Agent plan:
                1. Use query_openbis_objects using type = "SUBSTANCE" and filters = {{"empa_number": 704, "batch": 'a'}}
                2. Use query_openbis_experiments to list all experiments.
                3. Use query_openbis_objects to get all the data of the objects connected to the experiments
                4. Keep using query_openbis_objects for finding data about objects connected to other objects
                5. Filter for experiments containing that substance (like by the name or permId)
                6. Summarize relevant experiments.

            Tools used:
            - query_openbis_objects
            - query_openbis_experiments

            Today is {get_current_time()}.
        """

        self._tools = [
            query_openbis_experiments,
            # get_simulations,
            query_openbis_objects,
        ]

        llm_with_tools = self.llm_model.bind_tools(self._tools)

        def chatbot(state: State):
            # Inject the system prompt at the start (only once)
            messages = state["messages"]

            # Check if the system prompt is already there to avoid duplication
            if not any(
                isinstance(msg, SystemMessage) and msg.content == self.system_prompt
                for msg in messages
            ):
                messages = [SystemMessage(content=self.system_prompt)] + messages

            return {"messages": [llm_with_tools.invoke(messages)]}

        # Define the state machine
        # Add the nodes
        graph_builder = StateGraph(State)
        graph_builder.add_node("chatbot", chatbot)
        tool_node = ToolNode(tools=self._tools)
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
