from langgraph.graph import StateGraph, END
from .state import AgentState
from .planner import planner_node
from .executor import executor_node
from .reflection import reflection_node
from .summary import summary_node
from .ask_user import ask_user_node
from .respond import respond_node

# 1. Initialize Graph
workflow = StateGraph(AgentState)

# 2. Add Nodes
workflow.add_node("planner", planner_node)
workflow.add_node("respond", respond_node)
workflow.add_node("ask_user_node", ask_user_node)
workflow.add_node("executor", executor_node)
workflow.add_node("reflection", reflection_node)
workflow.add_node("summary", summary_node)

# 3. Set Entry Point
workflow.set_entry_point("planner")


# 4. Define Routing Logic
def route_after_planning(state: AgentState):
    """
    THREE-WAY ROUTING after Planner
    """
    action = state.get("next_action", "ask_user")

    if action == "respond":
        return "respond"
    elif action == "execute":
        return "executor"  # Tool execution

    # Fallback
    return "ask_user_node"


def route_after_reflection(state: AgentState):
    """
    Routing after Reflection
    """
    # Default is "summary" only if action is explicitly summary or unknown
    action = state.get("next_action", "summary")

    if action == "replan":
        return "planner"  # Try again
    elif action == "ask_user":
        return "ask_user_node"  # Blocker
    # elif action == "respond":
    #     return "respond"  # Logical Dead End -> Direct explanation

    return "summary"  # Success default


# 5. Register Edges
workflow.add_conditional_edges("planner", route_after_planning)
workflow.add_edge("executor", "reflection")
workflow.add_conditional_edges("reflection", route_after_reflection)

# End Points
# workflow.add_edge("respond", END)
workflow.add_edge("summary", END)
workflow.add_edge("ask_user_node", END)

# 6. Compile
app = workflow.compile()
