from datetime import datetime
from typing import List, Optional, Dict, Any
from peewee import DoesNotExist, IntegrityError, fn

from ..model.agent_message_model import AgentMessage


class AgentMessageDAO:
    """
    Data Access Object (DAO) for managing AgentMessage database operations.
    Handles the storage and retrieval of granular conversation events
    (user inputs, tool calls, agent thoughts) within a session.
    """

    @staticmethod
    def create_message(
            message_id: str,
            session_id: str,
            message_type: str,
            message_data: Dict[str, Any]
    ) -> Optional[AgentMessage]:
        """
        Creates a new message record within a specific session.

        Args:
            message_id (str): Unique internal ID (e.g., UUID).
            session_id (str): The business ID of the parent session.
            message_type (str): Type of event (user_input, agent_thought, etc.).
            message_data (dict): The actual content payload (JSON).

        Returns:
            AgentMessage: The created message instance if successful.
            None: If creation fails (e.g., session not found).
        """
        try:
            now = datetime.utcnow()

            # Create Record
            # Note: We assign the string 'session_id' to the ForeignKey field 'session'.
            # Peewee resolves this to the related AgentSession model.
            message = AgentMessage.create(
                message_id=message_id,
                session=session_id,
                message_type=message_type,
                message_data=message_data,
                created_at=now
            )
            return message

        except IntegrityError as e:
            print(f"[Error] Integrity error creating message '{message_id}': {e}")
            return None
        except Exception as e:
            print(f"[Error] Unexpected error creating message: {e}")
            return None

    @staticmethod
    def get_message_by_id(message_id: str) -> Optional[AgentMessage]:
        """
        Retrieves a single message by its unique internal ID.

        Args:
            message_id (str): The unique message identifier.

        Returns:
            AgentMessage: The found object or None.
        """
        try:
            return AgentMessage.get(AgentMessage.message_id == message_id)
        except DoesNotExist:
            print(f"[Warning] Message '{message_id}' not found.")
            return None

    @staticmethod
    def get_messages_by_session(
            session_id: str,
            page: int = 1,
            page_size: int = 20,
            types: Optional[List[str]] = None
    ) -> List[AgentMessage]:
        """
        Retrieves paginated history for a session, ordered by newest first.
        Essential for frontend chat interfaces (load previous messages).

        Args:
            session_id (str): The session to fetch history for.
            page (int): Page number (1-based index).
            page_size (int): Number of messages per page.
            types (List[str], optional): Filter by specific message types.
                                        (e.g., pass ['user_input', 'agent_final']
                                        to hide internal thoughts/tool calls).

        Returns:
            List[AgentMessage]: List of message objects.
        """
        try:
            # 1. Base Query
            query = AgentMessage.select().where(AgentMessage.session == session_id)

            # 2. Apply Type Filter (if provided)
            if types:
                query = query.where(AgentMessage.message_type.in_(types))

            # 3. Order By Newest First (Desc)
            # This ensures the frontend gets the most recent context first.
            query = query.order_by(AgentMessage.created_at.desc())

            # 4. Apply Pagination
            # Peewee's paginate method uses (page_number, items_per_page)
            paginated_query = query.paginate(page, page_size)

            return list(paginated_query)

        except Exception as e:
            print(f"[Error] Failed to retrieve messages for session '{session_id}': {e}")
            return []

    @staticmethod
    def get_messages_page(
            session_id: str,
            page: int = 1,
            page_size: int = 20,
            ascending: bool = False,
            types: Optional[List[str]] = None
    ) -> List[Dict[str, Any]]:
        """
        Retrieves a paginated slice of messages for frontend display.
        frontend needs pagination (page=1,2,3...) and may require both ascending or descending order.

        Args:
            session_id (str): The session identifier.
            page (int): Page number starting from 1.
            page_size (int): Number of messages per page.
            ascending (bool): Whether to return messages in chronological order.
            types (List[str], optional): Filter by specific message types.

        Returns:
            List[Dict]: List of message dictionaries for frontend rendering.
        """
        try:
            # 1. Compute pagination offsets
            offset = (page - 1) * page_size

            # 2. Determine ordering direction
            order_clause = (AgentMessage.created_at.asc()
                            if ascending
                            else AgentMessage.created_at.desc())

            # 3. Query the paginated messages
            query = (AgentMessage
                     .select()
                     .where(AgentMessage.session == session_id))

            if types:
                query = query.where(AgentMessage.message_type.in_(types))

            query = (query
                     .order_by(order_clause)
                     .offset(offset)
                     .limit(page_size))

            # 4. Convert ORM objects into frontend-friendly dicts
            results = []
            for msg in query:
                results.append({
                    "id": msg.id,
                    "type": msg.message_type,
                    "content": msg.message_data,
                    "timestamp": msg.created_at.isoformat()
                })

            return results

        except Exception as e:
            print(f"[Error] Failed to fetch paginated messages for '{session_id}': {e}")
            return []

    @staticmethod
    def get_agent_context_window(
            session_id: str,
            limit: int = 10
    ) -> List[Dict[str, Any]]:
        """
        Retrieves recent history formatted for the Agent's prompt context.
        Unlike the frontend view, this usually requires chronological order (Oldest -> Newest)
        so the LLM can read the conversation flow naturally.

        Args:
            session_id (str): The session identifier.
            limit (int): How many recent messages to include in the context window.

        Returns:
            List[Dict]: List of dictionaries ready for LLM injection.
        """
        try:
            # 1. Get the N most recent messages
            subquery = (AgentMessage
                        .select()
                        .where(AgentMessage.session == session_id)
                        .order_by(AgentMessage.created_at.desc())
                        .limit(limit))

            # 2. Re-order them to Chronological (Ascending) for the LLM
            # We fetch the subquery and reverse it using python or a surrounding query
            messages = sorted(subquery, key=lambda x: x.created_at)

            # 3. Construct Context Dicts
            context_list = []
            for msg in messages:
                context_list.append({
                    "role": msg.message_type,  # Logic might map 'user_input' -> 'user'
                    "content": msg.message_data,
                    "timestamp": msg.created_at.isoformat()
                })

            return context_list

        except Exception as e:
            print(f"[Error] Failed to build agent context for '{session_id}': {e}")
            return []

    @staticmethod
    def count_messages(session_id: str, types: Optional[List[str]] = None) -> int:
        """
        Counts total messages in a session. Useful for statistics.
        """
        try:
            query = AgentMessage.select().where(AgentMessage.session == session_id)
            if types:
                query = query.where(AgentMessage.message_type.in_(types))
            return query.count()
        except Exception:
            return 0
