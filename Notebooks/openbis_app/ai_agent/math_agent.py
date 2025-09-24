from agent_state import AgentState

# --- Math Agent ---
def math_agent(state: AgentState) -> str:
    """
    A math-solving agent that uses the LLM to process and solve math problems.

    Args:
        state (AgentState): Contains the user's query.

    Returns:
        dict: Updated state with the computed answer from the LLM.
    """
    print("--- Math Node ---")
    prompt = f"Solve this math problem and return only the answer: {state['user_query']}"
    response = llm.invoke(prompt)
    state['answer'] = response.content.strip()
    return state