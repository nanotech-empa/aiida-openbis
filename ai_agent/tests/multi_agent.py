import tools
from langchain_google_genai import ChatGoogleGenerativeAI
from zoneinfo import ZoneInfo
from datetime import datetime
from langgraph.checkpoint.memory import InMemorySaver
from typing import Annotated
from langgraph.graph import StateGraph, START, END, MessagesState
from langgraph.prebuilt import create_react_agent, InjectedState
from langchain_core.tools import tool, InjectedToolCallId
from langgraph.types import Command
from uuid import uuid4
import json
from typing import Optional
import warnings

warnings.filterwarnings("ignore")


def read_json(filename):
    with open(filename, "r") as file:
        return json.load(file)


def read_text_file(file_path: str) -> str:
    with open(file_path, "r") as file:
        return file.read()


def get_current_time() -> str:
    """
    Get the current time formatted nicely for display.
    """
    # Get current datetime
    now = datetime.now(ZoneInfo("Europe/Zurich"))
    formatted = now.strftime("%A, %B %d, %Y at %H:%M")
    return f"It is currently {formatted}."


def create_handoff_tool(*, agent_name: str, description: Optional[str] = None):
    name = f"transfer_to_{agent_name}"
    description = description or f"Ask {agent_name} for help."

    @tool(name, description=description)
    def handoff_tool(
        state: Annotated[MessagesState, InjectedState],
        tool_call_id: Annotated[str, InjectedToolCallId],
    ) -> Command:
        tool_message = {
            "role": "tool",
            "content": f"Successfully transferred to {agent_name}",
            "name": name,
            "tool_call_id": tool_call_id,
        }
        return Command(
            goto=agent_name,
            update={**state, "messages": state["messages"] + [tool_message]},
            graph=Command.PARENT,
        )

    return handoff_tool


LLM_API_KEY = read_json("/home/jovyan/api_keys/gemini_api.json")["api_key"]
LLM_MODEL = ChatGoogleGenerativeAI(
    model="models/gemini-2.5-flash", google_api_key=LLM_API_KEY
)


class OpenBISAgent:
    """
    An agent that can answer questions about openBIS inventories.
    It uses a LangChain model and tools to interact with the openBIS session.
    """

    def __init__(self):
        # Define the sub-agents
        general_agent_prompt = read_text_file("data/system_prompt.txt")
        general_agent = create_react_agent(
            model=LLM_MODEL,
            tools=[
                tools.get_openbis_objects,
                tools.get_openbis_object_by_permId,
                tools.get_openbis_objects_by_name,
                tools.get_openbis_objects_by_date,
            ],
            prompt=f"""
                {general_agent_prompt}
                INSTRUCTIONS:
                - Assist with general openBIS-related tasks, e.g., getting objects by type,
                permId, name, or date.
                - After you're done with your tasks, respond to the supervisor directly.
                - Respond ONLY with the results of your work, do NOT include ANY other text.
            """,
            name="general_agent",
        )

        specific_agent = create_react_agent(
            model=LLM_MODEL,
            tools=[
                tools.get_live_samples_by_attributes,
                tools.get_substances_by_attributes,
                # tools.get_crystals_by_attributes,
                # tools.get_2d_materials_by_attributes
            ],
            prompt="""
                You are an openBIS agent able to answer specific questions about openBIS inventories, e.g.,
                when the user wants to find substances by attributes, or get live samples.
                INSTRUCTIONS:
                - Assist ONLY with specific openBIS-related tasks.
                - After you're done with your tasks, respond to the supervisor directly.
                - Respond ONLY with the results of your work, do NOT include ANY other text.
            """,
            name="specific_agent",
        )

        # Handoffs
        assign_to_general_agent = create_handoff_tool(
            agent_name="general_agent",
            description="Assign task to a general agent.",
        )

        assign_to_specific_agent = create_handoff_tool(
            agent_name="specific_agent",
            description="Assign task to a specific agent.",
        )

        supervisor_agent = create_react_agent(
            model=LLM_MODEL,
            tools=[assign_to_general_agent, assign_to_specific_agent],
            prompt="""
                You are a supervisor managing two agents:
                - a general agent. Assign general openBIS tasks to this agent, e.g. getting objects by type, permId,
                name, date, or description.
                - a specific agent. Assign specific openBIS tasks to this agent, e.g. getting substances by attributes,
                getting crystals by attributes, getting 2D materials by attributes, or getting live samples.
                After any agent responds, you must always generate a final user-facing message.
                Never leave your output empty.
            """,
            name="supervisor",
        )

        # Define the state machine
        # Add the nodes
        memory = InMemorySaver()
        self.agent = (
            StateGraph(MessagesState)
            .add_node(
                supervisor_agent, destinations=("general_agent", "specific_agent", END)
            )
            .add_node(general_agent)
            .add_node(specific_agent)
            .add_edge(START, "supervisor")
            .add_edge("general_agent", "supervisor")
            .add_edge("specific_agent", "supervisor")
            .compile(checkpointer=memory)
        )

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
            {"messages": [{"role": "user", "content": user_prompt}]},
            self.agent_config,
            stream_mode="values",
        )
        response = []
        for event in events:
            # print(event["messages"][-1])
            # pretty_print_messages(event, last_message=True)
            response.append(event["messages"][-1])
        return response

    def get_graph_replay(self):
        to_replay = None
        for state in self.agent.get_state_history(self.agent_config):
            print("Num Messages: ", len(state.values["messages"]), "Next: ", state.next)
            print("-" * 80)
            if len(state.values["messages"]) == 10:
                # We are somewhat arbitrarily selecting a specific state based on the number of chat messages in the state.
                to_replay = state

        return to_replay
