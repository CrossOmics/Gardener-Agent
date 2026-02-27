"""
Agent Message table for storing conversation events.

Represents individual messages within an agent session, including user inputs,
agent reasoning steps, tool calls, tool results, and final responses.

Messages are linked to agent_sessions via session_id and stored as structured
JSON payloads to support flexible, multi-step ReAct-style interactions.
"""
from peewee import (
    AutoField,
    CharField,
    DateTimeField,
    ForeignKeyField
)
from datetime import datetime

from playhouse.sqlite_ext import JSONField

from .base_model import BaseModel
from .agent_session_model import AgentSession


class AgentMessage(BaseModel):
    """
    Agent Message table for storing conversation events.

    Each record represents a single message or internal step (user input,
    agent thought, tool call/result, or final answer) within an agent session.
    """

    # Primary Key
    id = AutoField(
        primary_key=True,
        help_text="Auto-increment primary key"
    )

    # Internal Message ID
    message_id = CharField(
        max_length=128,
        unique=True,
        index=True,
        null=False,
        help_text="Unique internal message identifier"
    )

    # Session Foreign Key (business key)
    session = ForeignKeyField(
        AgentSession,
        field='session_id',
        column_name='session_id',

        on_delete='CASCADE',
        backref='message',
        help_text="Reference to the parent project's business ID"
    )

    # Message Type
    message_type = CharField(
        max_length=32,
        null=False,
        help_text=(
            "Message type: user_input, agent_thought, agent_plan, "
            "agent_action, agent_observation, tool_call, tool_result, agent_final"
        )
    )

    # Message Payload
    message_data = JSONField(
        null=False,
        help_text="Message content stored as structured JSON"
    )

    # Temporal Field
    created_at = DateTimeField(
        default=datetime.utcnow,
        null=False,
        help_text="Message creation timestamp for ordering"
    )

    class Meta:
        table_name = "agent_messages"
        indexes = (
            # Optimized for fetching conversation history
            (("session", "created_at"), False),
            (("session", "message_type"), False),
        )

    def __repr__(self) -> str:
        return (
            f"<AgentMessage id={self.id} "
            f"message_id='{self.message_id}' "
            f"type='{self.message_type}'>"
        )
