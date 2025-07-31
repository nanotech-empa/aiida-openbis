from ai_agent import openbis_utils
from ai_agent.tools import inventory_tools

import asyncio
import ipywidgets as ipw
from IPython.display import display
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
from langchain_core.messages import SystemMessage
from uuid import uuid4

def get_current_time() -> str:
    """
    Get the current time formatted nicely for display.
    """
    # Get current datetime
    now = datetime.now(ZoneInfo("Europe/Zurich"))
    formatted = now.strftime("%A, %B %d, %Y at %H:%M")
    return f"It is currently {formatted}."

class State(TypedDict):
    messages: Annotated[list, add_messages]

class OpenBISAgent():
    """
    An agent that can answer questions about openBIS inventories.
    It uses a LangChain model and tools to interact with the openBIS session.
    """
    def __init__(self, google_api_key):
        self.google_api_key = google_api_key
        self.llm_model = ChatGoogleGenerativeAI(model = "models/gemini-2.5-flash-lite", google_api_key = self.google_api_key)
        
        self.system_prompt = f"""
            You are a helpful assistant that can answer questions about experiments, simulations, 
            and all the inventories inside openBIS.
            
            In our openBIS instance there are multiple types, such as:
            SUBSTANCE: substances (also known as precursors) and chemicals. However there is a difference between them
            precursors are located in collection /MATERIALS/MOLECULES/PRECURSOR_COLLECTION and chemicals are located
            in collection /MATERIALS/RAW_MATERIALS/CHEMICAL_COLLECTION
            MOLECULE: molecule that is the basis for the different substances and chemicals
            PUBLICATION: publications
            GRANT: grants
            ORGANISATION: companies, organisations, and institutions
            PERSON: people
            LOCATION: locations such as rooms
            INSTRUMENT: instruments found in the labs
            INSTRUMENT.STM: instruments found in the labs specialised for STM measurements
            ANALYSER: analyser component
            CHAMBER: chamber component
            COMPONENT: general component
            ELECTRONICS: electronics component
            EVAPORATOR: evaporator component
            EVAPORATOR_SLOT: evaporator slot component which is normally part of the evaporator
            ION_GAUGE: ion gauge component
            ION_PUMP: ion pump component
            PBN_STAGE: PBN stage component
            SCROLL_PUMP: scroll pump component
            SPUTTER_GUN: sputter gun component
            STM_AFM_TIP: stm/afm tip component
            THYRACONT: thyracont component
            TURBO_PUMP: turbo pump component
            VALVE: valve component
            2D_LAYER_MATERIAL: 2D layer material
            WAFER_SAMPLE: wafer sample
            WAFER: wafer
            WIRE: wire
            SAMPLE: samples prepared in the labs using annealing, deposition, sputtering tasks, etc.
            ATOMISTIC_MODEL: atomistic models developed in the computational simulations
            AIIDA_NODE: AiiDA nodes representing the nodes inside AiiDAlab 
            CRYSTAL: crystals (also known as slabs)
            CRYSTAL_CONCEPT: theorical concept of the crystals used in the labs
            SOFTWARE: software used in the analysis of measurements
            CODE: scripts/codes used in the analysis of measurements or for performing simulations
            
            Today is {get_current_time()}.
        """
        
        self._tools = [
            # General openBIS objects tools
            inventory_tools.get_openbis_objects,
            inventory_tools.get_openbis_object_by_permId,
            inventory_tools.get_openbis_objects_by_name,
            inventory_tools.get_openbis_objects_by_date,
            
            # Substances tools
            inventory_tools.get_substances_by_attributes,
            inventory_tools.get_crystals_by_attributes,
            inventory_tools.get_2d_materials_by_attributes,
        ]
        
        llm_with_tools = self.llm_model.bind_tools(self._tools)
        def chatbot(state: State):
            # Inject the system prompt at the start (only once)
            messages = state["messages"]
            
            # Check if the system prompt is already there to avoid duplication
            if not any(isinstance(msg, SystemMessage) and msg.content == self.system_prompt for msg in messages):
                messages = [SystemMessage(content=self.system_prompt)] + messages

            return {"messages": [llm_with_tools.invoke(messages)]}
        
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
