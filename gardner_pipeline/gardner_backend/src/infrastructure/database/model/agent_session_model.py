from peewee import (
    AutoField,
    CharField,
    DateTimeField,
    TextField, ForeignKeyField
)
from datetime import datetime

from .project_meta_model import ProjectMeta
from .base_model import BaseModel


class AgentSession(BaseModel):
    """
    Agent Session table for managing conversational contexts.

    Each session represents a single logical conversation between the user
    and an agent, grouping all related messages, tool calls, and reasoning
    steps under a shared session_id.
    """

    # Primary Key
    id = AutoField(
        primary_key=True,
        help_text="Auto-increment primary key"
    )

    # Business Session ID
    session_id = CharField(
        max_length=128,
        unique=True,
        index=True,
        null=False,
        help_text="Conversation session identifier (e.g., session_proj_001_analysis)"
    )

    # Foreign Key links to ProjectMeta's project id
    project = ForeignKeyField(
        ProjectMeta,
        field='project_id',
        column_name='project_id',

        on_delete='CASCADE',
        backref='session',
        help_text="Reference to the parent project's business ID"
    )
    # Agent Identity
    session_name = CharField(
        max_length=256,
        null=False,
        help_text="Agent session name"
    )
    # Agent Identity
    agent_name = CharField(
        max_length=64,
        null=False,
        help_text="Agent model name (e.g., gpt-4.1, gemini-pro)"
    )

    # Session Lifecycle Status
    status = CharField(
        max_length=32,
        null=False,
        default="INIT",
        help_text=(
            "Session execution status: "
            "INIT, RUNNING, WAITING, COMPLETED, FAILED, CANCELLED"
        )
    )

    # Temporal Fields
    created_at = DateTimeField(
        default=datetime.utcnow,
        null=False,
        help_text="Session creation timestamp"
    )

    updated_at = DateTimeField(
        default=datetime.utcnow,
        null=False,
        help_text="Last status update timestamp"
    )

    class Meta:
        table_name = "agent_sessions"
        indexes = (
            # Common query: all sessions under a project
            (("project_id", "created_at"), False),
        )

    def __repr__(self) -> str:
        return (
            f"<AgentSession id={self.id} "
            f"session_id='{self.session_id}' "
            f"status='{self.status}'>"
        )
