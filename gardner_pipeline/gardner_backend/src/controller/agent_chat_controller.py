from typing import List, Dict, Any
from fastapi import APIRouter, Depends, HTTPException, status, Query, Path

from core.task_executor import task_executor
from dto.request.agent_message_request import AgentMessageRequest
from dto.request.create_session_request import CreateSessionRequest
from dto.request.rename_session_request import RenameSessionRequest
from service.agent_chat_service import AgentChatService

router = APIRouter(
    prefix="/api/v1/agent/history",
    tags=["Biological Agent Chat History"]
)


# Session Endpoints
@router.post("/sessions/create", status_code=status.HTTP_201_CREATED)
async def create_chat_session(
        request: CreateSessionRequest,
        service: AgentChatService = Depends()
):
    def _task():
        return service.start_new_session(
            project_id=request.project_id,
            agent_name=request.agent_name,
            initial_name=request.initial_name
        )

    try:
        return await task_executor.run_in_thread(_task)
    except ValueError as e:
        raise HTTPException(status_code=404, detail=str(e))
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/project/session", response_model=List[Dict[str, Any]])
async def list_project_sessions(
        project_id: str = Query(...),
        service: AgentChatService = Depends()
):
    """
    Retrieves all conversation records for a project, ordered from newest to oldest.
    """
    def _task():
        return service.list_project_sessions(project_id)

    return await task_executor.run_in_thread(_task)


@router.patch("/sessions/rename/{session_id}")
async def rename_chat_session(
        request: RenameSessionRequest,
        session_id: str = Path(...),
        service: AgentChatService = Depends()
):
    def _task():
        if not service.rename_session(session_id, request.new_name):
            raise HTTPException(status_code=404, detail="Session not found")
        return {"msg": "Successfully Renamed"}

    return await task_executor.run_in_thread(_task)


@router.delete("/sessions/delete/{session_id}", status_code=status.HTTP_204_NO_CONTENT)
async def delete_chat_session(
        session_id: str = Path(...),
        service: AgentChatService = Depends()
):
    def _task():
        if not service.delete_session(session_id):
            raise HTTPException(status_code=404, detail="Session not found")

    await task_executor.run_in_thread(_task)


# Unified Logging Endpoint
@router.post("/{session_id}/message/create", status_code=status.HTTP_201_CREATED)
async def create_message(
        request: AgentMessageRequest,
        session_id: str = Path(..., description="Session Business ID"),
        service: AgentChatService = Depends()
):
    """
    Unified endpoint to log User inputs, Agent thoughts, Tool calls, etc.
    Dispatches to the correct service method based on `message_type`.
    """
    msg_type = request.message_type
    payload = request.data

    def _task():
        # Dispatcher Logic: Maps generic JSON to specific Service methods
        if msg_type == "user_input":
            msg_id = service.log_user_input(
                session_id,
                payload.get("text", "")
            )
        elif msg_type == "agent_thought":
            msg_id = service.log_agent_thought(
                session_id,
                payload.get("thought", "")
            )
        elif msg_type == "tool_call":
            msg_id = service.log_tool_call(
                session_id,
                payload.get("tool_name", "unknown"),
                payload.get("arguments", {})
            )
        elif msg_type == "tool_result":
            msg_id = service.log_tool_result(
                session_id,
                payload.get("tool_name", "unknown"),
                payload.get("result_data"),
                payload.get("is_error", False)
            )
        elif msg_type == "agent_final":
            msg_id = service.log_agent_response(
                session_id,
                payload.get("text", ""),
                payload.get("attachments", [])
            )
        else:
            msg_id = service.log_event(session_id, msg_type, payload)

        return {"msg_id": msg_id, "status": "logged", "type": msg_type}

    try:
        return await task_executor.run_in_thread(_task)
    except KeyError as e:
        raise HTTPException(status_code=400, detail=f"Missing field in payload: {e}")
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Logging failed: {str(e)}")


# History Retrieval
@router.get("/{session_id}/history/all")
async def get_session_history(
        session_id: str = Path(...),
        page: int = Query(1, ge=1),
        size: int = Query(20, ge=1, le=100),
        service: AgentChatService = Depends()
):
    """
    Get paginated history for UI rendering.
    """
    def _task():
        return service.get_session_history_ui(session_id, page, size)

    return await task_executor.run_in_thread(_task)


@router.get("/{session_id}/interaction")
async def get_chat_only_history(
        session_id: str = Path(...),
        page: int = Query(1, ge=1),
        size: int = Query(20, ge=1, le=100),
        service: AgentChatService = Depends()
):
    """
    Get paginated history containing only user inputs and agent final responses.
    """
    def _task():
        return service.get_chat_only_history(session_id, page, size)

    return await task_executor.run_in_thread(_task)
