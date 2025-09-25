import tools
from langchain_google_genai import ChatGoogleGenerativeAI
from zoneinfo import ZoneInfo
from datetime import datetime
from langgraph.checkpoint.memory import InMemorySaver
from langgraph.prebuilt import ToolNode, tools_condition
from typing import Annotated
from typing_extensions import TypedDict
from langgraph.graph import StateGraph, START, END
from langgraph.graph.message import add_messages
from langchain_core.messages import SystemMessage
from langgraph.prebuilt import create_react_agent
from uuid import uuid4
import json
from typing import Literal
import warnings
warnings.filterwarnings("ignore")

def read_json(filename):
    with open(filename, 'r') as file:
        return json.load(file)

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

class AgentState(TypedDict):
    user_prompt: str
    answer: str
    messages: Annotated[list, add_messages]

LLM_API_KEY = read_json("/home/jovyan/api_keys/gemini_api.json")["api_key"]
LLM_MODEL = ChatGoogleGenerativeAI(model = "models/gemini-2.5-flash", google_api_key = LLM_API_KEY)

def general_tools_agent(state: AgentState) -> dict:
    """
    Executes a ReAct-style agent that processes a user query.
    This agent is used for answering general queries about openBIS objects, e.g.
    getting objects of a certain type, or getting an object by its permId, name, date, or description.

    This function takes the current state (which includes the user's question)
    creates an agent using the LLM and a set of predefined tools, then runs the agent to get a response.
    """
    agent_system_prompt = read_text_file("data/system_prompt.txt")
    agent = create_react_agent(
        model = LLM_MODEL,
        tools = [
            tools.get_openbis_objects,
            tools.get_openbis_object_by_permId,
            tools.get_openbis_objects_by_name,
            tools.get_openbis_objects_by_date
        ],
        prompt = agent_system_prompt
    )
    result = agent.invoke({"messages": state["user_prompt"]})
    print(result["messages"][-1].content)
    return {"answer": result["messages"][-1].content}

def specific_tools_agent(state: AgentState) -> dict:
    """
    Executes a ReAct-style agent that processes a user query.
    This agent is used for answering specific queries about openBIS objects, e.g.
    getting substances by attributes, getting crystals by attributes, 
    getting 2D materials by attributes, or getting live samples.

    This function takes the current state (which includes the user's question)
    creates an agent using the LLM and a set of predefined tools, then runs the agent to get a response.
    """
    agent = create_react_agent(
        model = LLM_MODEL,
        tools = [
            tools.get_live_samples_by_attributes,
            tools.get_substances_by_attributes,
            # tools.get_crystals_by_attributes,
            # tools.get_2d_materials_by_attributes
        ]
    )
    result = agent.invoke({"messages": state["user_prompt"]})
    return {"answer": result["messages"][-1].content}

agent_docs = {
    "general_tools_agent": general_tools_agent.__doc__,
    "specific_tools_agent": specific_tools_agent.__doc__
}

def supervisor_agent(state: AgentState) -> dict:
    """
    Supervisor agent that decides whether to answer directly or ask sub-agents.
    """
    query = state["user_prompt"]

    decision_prompt = f"""
    You are the supervisor agent. The user asked: "{query}".
    
    You have two options:
    1. Answer directly if it is general knowledge or does not require openBIS tools.
    2. If the question is about openBIS data, choose which sub-agent should handle it:
       - general_tools_agent: {agent_docs['general_tools_agent']}
       - specific_tools_agent: {agent_docs['specific_tools_agent']}
    
    Respond with one of:
    - "direct:<your answer here>"
    - "general:<reason why general agent should be used>"
    - "specific:<reason why specific agent should be used>"
    """
    decision = LLM_MODEL.invoke(decision_prompt).content.strip()

    if decision.startswith("direct:"):
        answer = decision.replace("direct:", "").strip()
        return {"answer": answer, "messages": [{"role": "assistant", "content": answer}]}

    elif decision.startswith("general:"):
        sub_result = general_tools_agent(state)
        final_answer = f"I asked the general agent: {sub_result['answer']}"
        return {"answer": final_answer, "messages": [{"role": "assistant", "content": final_answer}]}

    elif decision.startswith("specific:"):
        sub_result = specific_tools_agent(state)
        final_answer = f"I asked the specific agent: {sub_result['answer']}"
        return {"answer": final_answer, "messages": [{"role": "assistant", "content": final_answer}]}

    else:
        return {"answer": "I am not sure how to handle this query.", 
                "messages": [{"role": "assistant", "content": "I am not sure how to handle this query."}]}

class OpenBISAgent():
    """
    An agent that can answer questions about openBIS inventories.
    It uses a LangChain model and tools to interact with the openBIS session.
    """
    def __init__(self):
        # Define the state machine
        # Add the nodes
        workflow = StateGraph(AgentState)
        workflow.add_node("supervisor_agent", supervisor_agent)
        workflow.add_edge(START, "supervisor_agent")
        workflow.add_edge("supervisor_agent", END)
        
        memory = InMemorySaver()
        self.agent = workflow.compile(checkpointer=memory)
        self.agent_config = {"configurable": {"thread_id": uuid4()}}
    
    def ask_question(self, user_prompt: str):
        """
        Ask a question to the agent and get a response.
        
        Args:
            user_prompt (str): The question to ask the agent.
        
        Returns:
            str: The response from the agent.
        """
        events = self.agent.stream(
            {"user_prompt": user_prompt, "messages": [{"role": "user", "content": user_prompt}]},
            self.agent_config,
            stream_mode="values",
        )
        response = []
        for event in events:
            response.append(event["messages"][-1])
        return response