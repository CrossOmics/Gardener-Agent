from typing import List, Dict, Any
from fastapi import Depends

from infrastructure.database.dao.agent_session_dao import AgentSessionDAO
from infrastructure.database.dao.agent_message_dao import AgentMessageDAO
from util.id_generate_utils import generate_business_id, generate_session_id_from_project


class AgentChatService:
    """
    Service layer for managing biological agent conversations.

    Responsibilities:
    1. Session Management: Create, list, rename, and delete chat sessions.
    2. Message Tracing: Log distinct steps of the agent loop (User Input -> Thought -> Tool -> Result -> Final).
    3. Context Retrieval: Fetch history for both the UI (display) and the LLM (context window).
    """

    def __init__(
            self,
            session_dao: AgentSessionDAO = Depends(),
            message_dao: AgentMessageDAO = Depends()
    ):
        self.session_dao = session_dao
        self.message_dao = message_dao

    # 1. Session Lifecycle Management
    def start_new_session(self, project_id: str, agent_name: str = "BioAgent-Pro",
                          initial_name: str = "New Analysis") -> Dict[str, Any]:
        """
        Initializes a new chat session linked to a specific bio-project.
        """
        # Generate a unique business ID for the session based on project id
        session_id = generate_session_id_from_project(project_id, "session")

        session = self.session_dao.create_session(
            session_id=session_id,
            project_id=project_id,
            agent_name=agent_name,
            session_name=initial_name,
            status="INIT"
        )

        if not session:
            raise ValueError(f"Failed to create session for project {project_id}")

        return {
            "session_id": session.session_id,
            # Peewee internal field access
            "project_id": session.project,
            "name": session.session_name,
            "status": session.status,
            "created_at": session.created_at
        }

    def list_project_sessions(self, project_id: str) -> List[Dict[str, Any]]:
        """
        Lists all conversation history for a project, useful for a sidebar menu.
        """
        sessions = self.session_dao.get_sessions_by_project(project_id)
        return [
            {
                "session_id": s.session_id,
                "name": s.session_name,
                "agent": s.agent_name,
                "updated_at": s.updated_at
            }
            for s in sessions
        ]

    def rename_session(self, session_id: str, new_name: str) -> bool:
        """
        Renames a session (e.g., after the agent summarizes the topic).
        """
        return self.session_dao.update_session_name(session_id, new_name)

    def delete_session(self, session_id: str) -> bool:
        """
        Permanently deletes a session and its history.
        """
        return self.session_dao.delete_session(session_id)

    # 2. Message Logging (The Agent "Trace")
    def log_event(self, session_id: str, message_type: str, payload: Dict[str, Any]) -> str:
        """
        Internal helper to generate an ID and save a message record.
        """
        message_id = generate_business_id(f"m_{message_type}")

        msg = self.message_dao.create_message(
            message_id=message_id,
            session_id=session_id,
            message_type=message_type,
            message_data=payload
        )

        if not msg:
            raise RuntimeError(f"Failed to log event type {message_type}")

        # Update session status to indicate activity
        self.session_dao.update_session_status(session_id, "RUNNING")

        return message_id

    def log_user_input(self, session_id: str, text: str) -> str:
        """
        Logs the user's initial instruction or follow-up question.
        """
        return self.log_event(session_id, "user_input", {"content": text})

    def log_agent_thought(self, session_id: str, thought_text: str) -> str:
        """
        Logs the agent's internal reasoning (CoT) before acting.
        Frontend usually renders this in a collapsible 'Thinking...' box.
        """
        return self.log_event(session_id, "agent_thought", {"thought": thought_text})

    def log_tool_call(self, session_id: str, tool_name: str, arguments: Dict[str, Any]) -> str:
        """
        Logs the agent's decision to call an external tool (e.g., RunClustering).
        """
        payload = {
            "tool": tool_name,
            "arguments": arguments,
            "status": "initiated"
        }
        return self.log_event(session_id, "tool_call", payload)

    def log_tool_result(self, session_id: str, tool_name: str, result_data: Any, is_error: bool = False) -> str:
        """
        Logs the output returned by the tool.
        Result data can be a simple string or complex JSON (e.g., file paths, stats).
        """
        payload = {
            "tool": tool_name,
            "output": result_data,
            "status": "error" if is_error else "success"
        }
        return self.log_event(session_id, "tool_result", payload)

    def log_agent_response(self, session_id: str, text: str, attachments: List[str] = None) -> str:
        """
        Logs the final answer provided to the user.
        """
        payload = {
            "content": text,
            "attachments": attachments or []  # List of file paths or image URLs
        }
        # Mark session as waiting for next user input
        self.session_dao.update_session_status(session_id, "WAITING")
        return self.log_event(session_id, "agent_final", payload)

    # Context & History Retrieval
    def get_session_history_ui(self, session_id: str, page: int = 1, size: int = 50) -> Dict[str, Any]:
        """
        Fetches paginated history formatted for the Frontend UI.
        Standardizes the output structure.
        """
        # Fetch raw records
        raw_msgs = self.message_dao.get_messages_page(
            session_id=session_id,
            page=page,
            page_size=size,
            ascending=False  # UI usually loads latest first
        )

        # Get total count for pagination metadata
        total = self.message_dao.count_messages(session_id)

        return {
            "session_id": session_id,
            "total_messages": total,
            "page": page,
            "size": size,
            "messages": raw_msgs  # List of {id, role, content, timestamp}
        }

    def get_chat_only_history(self, session_id: str, page: int = 1, size: int = 50) -> Dict[str, Any]:
        """
        Fetches paginated history containing only user inputs and agent final responses.
        """
        target_types = ["user_input", "agent_final"]
        
        # Fetch raw records
        raw_msgs = self.message_dao.get_messages_page(
            session_id=session_id,
            page=page,
            page_size=size,
            ascending=False,
            types=target_types
        )

        # Get total count for pagination metadata
        total = self.message_dao.count_messages(session_id, types=target_types)

        return {
            "session_id": session_id,
            "total_messages": total,
            "page": page,
            "size": size,
            "messages": raw_msgs
        }

    def get_agent_memory_window(self, session_id: str, token_limit_heuristic: int = 10) -> List[Dict[str, Any]]:
        """
        Fetches the recent conversation turn-by-turn for the LLM context window.
        Uses 'get_agent_context_window' from DAO which sorts chronologically.
        """
        return self.message_dao.get_agent_context_window(session_id, limit=token_limit_heuristic)
