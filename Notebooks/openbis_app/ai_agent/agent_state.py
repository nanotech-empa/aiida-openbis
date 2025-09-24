from typing import TypedDict 

# --- Define State ---
class AgentState(TypedDict):
    user_query: str
    answer: str