from datetime import datetime
from typing import List, Optional
from peewee import DoesNotExist, IntegrityError

from ..model.agent_session_model import AgentSession
from ..model.project_meta_model import ProjectMeta


class AgentSessionDAO:
    """
    Data Access Object (DAO) for managing AgentSession database operations.
    Handles the lifecycle of conversational contexts (sessions) associated with
    biological analysis projects.
    """

    @staticmethod
    def create_session(
            session_id: str,
            project_id: str,
            agent_name: str,
            session_name: Optional[str] = None,  # Added optional session_name
            status: str = "INIT"
    ) -> Optional[AgentSession]:
        """
        Creates a new conversation session record.

        Args:
            session_id (str): Unique business identifier (e.g., 'session_proj_001_analysis').
            project_id (str): The parent project's business ID to link this session to.
            agent_name (str): Name of the model/agent being used (e.g., 'gpt-4.1').
            session_name (str, optional): Human-readable name. Defaults to "Bio Assistant".
            status (str, optional): Initial status. Defaults to "INIT".

        Returns:
            AgentSession: The created session instance if successful.
            None: If creation fails (e.g., duplicate session_id or invalid project_id).
        """
        try:
            # 1. Verify project exists
            project_exists = ProjectMeta.select().where(ProjectMeta.project_id == project_id).exists()
            if not project_exists:
                print(f"[Error] Failed to create session: Project '{project_id}' does not exist.")
                return None

            now = datetime.utcnow()

            # Default naming logic: Use provided name or fallback to default
            final_name = session_name if session_name and session_name.strip() else "Bio Assistant"
            project_obj = ProjectMeta.get_or_none(ProjectMeta.project_id == project_id)
            # 2. Create Record
            session = AgentSession.create(
                session_id=session_id,
                project_id=project_obj,
                agent_name=agent_name,
                session_name=final_name,
                status=status,
                created_at=now,
                updated_at=now
            )
            return session

        except IntegrityError as e:
            print(f"[Error] Failed to create session '{session_id}': {e}")
            return None
        except Exception as e:
            print(f"[Error] Unexpected error creating agent session: {e}")
            return None

    @staticmethod
    def get_session_by_id(session_id: str) -> Optional[AgentSession]:
        """
        Retrieves a session by its unique business ID.

        Args:
            session_id (str): The unique string identifier.

        Returns:
            AgentSession: The found session instance.
            None: If not found.
        """
        try:
            return AgentSession.get(AgentSession.session_id == session_id)
        except DoesNotExist:
            print(f"[Warning] Session '{session_id}' not found.")
            return None

    @staticmethod
    def get_sessions_by_project(project_id: str) -> List[AgentSession]:
        """
        Retrieves all conversation sessions linked to a specific project.
        Ordered by creation time (newest first).

        Args:
            project_id (str): The project's business identifier.

        Returns:
            List[AgentSession]: List of session objects.
        """
        try:
            return list(
                AgentSession.select()
                .where(AgentSession.project == project_id)
                .order_by(AgentSession.created_at.desc())
            )
        except Exception as e:
            print(f"[Error] Failed to retrieve sessions for project '{project_id}': {e}")
            return []

    @staticmethod
    def update_session_status(session_id: str, new_status: str) -> bool:
        """
        Updates the execution status of a session.
        Automatically updates the 'updated_at' timestamp.

        Args:
            session_id (str): The unique session identifier.
            new_status (str): The new status code (e.g., 'RUNNING', 'COMPLETED').

        Returns:
            bool: True if update succeeded, False otherwise.
        """
        try:
            session = AgentSession.get(AgentSession.session_id == session_id)
            session.status = new_status
            session.updated_at = datetime.utcnow()
            session.save()
            return True
        except DoesNotExist:
            print(f"[Error] Cannot update status: Session '{session_id}' not found.")
            return False
        except Exception as e:
            print(f"[Error] Failed to update session '{session_id}': {e}")
            return False

    @staticmethod
    def update_session_name(session_id: str, new_name: str) -> bool:
        """
        Updates the human-readable name of a session.
        Useful when the user renames a chat or the agent auto-summarizes a title.

        Args:
            session_id (str): The unique session identifier.
            new_name (str): The new name for the session.

        Returns:
            bool: True if update succeeded.
        """
        try:
            if not new_name or not new_name.strip():
                print("[Warning] Cannot update session name to empty string.")
                return False

            session = AgentSession.get(AgentSession.session_id == session_id)
            session.session_name = new_name.strip()
            session.updated_at = datetime.utcnow()
            session.save()
            return True
        except DoesNotExist:
            print(f"[Error] Cannot rename: Session '{session_id}' not found.")
            return False
        except Exception as e:
            print(f"[Error] Failed to rename session '{session_id}': {e}")
            return False

    @staticmethod
    def delete_session(session_id: str) -> bool:
        """
        Deletes a session record.
        Note: Due to 'ON DELETE CASCADE' in ChatMessage model (if implemented),
        deleting a session will likely delete all associated messages.

        Args:
            session_id (str): The unique session identifier.

        Returns:
            bool: True if deletion was successful.
        """
        try:
            session = AgentSession.get(AgentSession.session_id == session_id)
            session.delete_instance()
            print(f"[Info] Session '{session_id}' deleted successfully.")
            return True
        except DoesNotExist:
            print(f"[Warning] Cannot delete: Session '{session_id}' does not exist.")
            return False
        except Exception as e:
            print(f"[Error] Failed to delete session '{session_id}': {e}")
            return False
