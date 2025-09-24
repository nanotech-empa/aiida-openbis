import openbis_utils
import inventory_tools

import asyncio
import ipywidgets as ipw
from IPython.display import display
from zoneinfo import ZoneInfo
from datetime import datetime
from functools import lru_cache
from langchain_core.tools import tool
from langgraph.checkpoint.memory import InMemorySaver
from langgraph.prebuilt import ToolNode, tools_condition, create_react_agent
from langchain_openai import ChatOpenAI
from langchain_google_genai import ChatGoogleGenerativeAI
from langchain_ollama import ChatOllama
from typing import List, Dict, Literal, Annotated
from typing_extensions import TypedDict
from langgraph.graph import StateGraph, START, END
from langgraph.graph.message import add_messages
from langchain_core.messages import SystemMessage
from uuid import uuid4

def read_text_file(file_path: str) -> str:
    with open(file_path, 'r') as file:
        return file.read()

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
    def __init__(self, llm_api_key):
        self.llm_api_key = llm_api_key
        # self.llm_model = ChatOllama(
        #     model = "qwen3:4b", 
        #     base_url = "http://host.docker.internal:11434", 
        #     reasoning = False,
        #     temperature = 0.7,
        #     top_p = 0.8,
        #     top_k = 20
        # ) # Too slow on my laptop
        self.llm_model = ChatGoogleGenerativeAI(model = "models/gemini-2.5-flash", google_api_key = self.llm_api_key)
        # self.llm_model = ChatOpenAI(
        #     model = "x-ai/grok-4-fast:free",
        #     base_url = "https://openrouter.ai/api/v1",
        #     openai_api_key = self.llm_api_key
        # )
        
        self.system_prompt = read_text_file("ai_agent/data/system_prompt.txt")

        self._tools = [
            # General openBIS objects tools
            inventory_tools.get_openbis_objects,
            inventory_tools.get_openbis_object_by_permId,
            inventory_tools.get_openbis_objects_by_name,
            inventory_tools.get_openbis_objects_by_date,
            
            # Specific tools
            inventory_tools.get_live_samples_by_attributes,
            inventory_tools.get_substances_by_attributes,
            # inventory_tools.get_crystals_by_attributes,
            # inventory_tools.get_2d_materials_by_attributes,
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
