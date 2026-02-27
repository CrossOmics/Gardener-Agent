import uvicorn
import sys
import io
from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
import shelve
from typing import Optional

from base.constants.local_config import CONFIG_DB, CUR_LLM
from base.dtos.request.agent_request import AgentRequest
from base.dtos.request.set_key_request import SetKeyRequest
from base.dtos.request.verify_key_request import VerifyKeyRequest
from base.dtos.response.agent_response import AgentResponse
from utils.api_utils import validate_api_key_online

# Force UTF-8 Encoding for Windows Console
if sys.stdout:
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
if sys.stderr:
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')

from orchestra.graph import app as agent_graph
from langchain_core.messages import HumanMessage
from loguru import logger

from base.base_api_manager import CredentialManager
from base.base_model_manager import reload_llm
from base.eums.model_type import APIKey

app = FastAPI(title="AI Agent Server")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# --- Initialize Credential Manager ---
credential_manager = CredentialManager(app_name="gardener_agent")


# --- Chat Endpoint ---
@app.post("/api/agent/chat", response_model=AgentResponse)
async def chat_endpoint(request: AgentRequest):
    """
    Main entry point for the Agent.
    Handles chat requests and executes the LangGraph workflow.
    """
    try:
        # 1. Prepare initial state for LangGraph
        initial_state = {
            "messages": [HumanMessage(content=request.message)],
            "project_id": request.project_id,
            "dataset_id": request.dataset_id,
            "snapshot_id": request.snapshot_id,
            "session_id": request.session_id,
            "iteration_count": 0,
            "execution_results": [],
            "filename": request.filename,
            "total_input_tokens": 0,
            "total_output_tokens": 0,
            "total_tool_calls": 0,
            "successful_tool_calls": 0,
        }

        # 2. Run the Workflow (This uses the LLM defined in config)
        # If API keys are missing, MockLLM will return a safe error message here.
        final_state = await agent_graph.ainvoke(initial_state)

        # 3. Extract the final response from the AI
        ai_message = final_state["messages"][-1]

        # Handle structured content (list of blocks)
        reply_content = ai_message.content
        if isinstance(reply_content, list):
            # Extract text from blocks like {'type': 'text', 'text': '...'}
            text_parts = []
            for block in reply_content:
                if isinstance(block, dict) and block.get("type") == "text":
                    text_parts.append(block.get("text", ""))
                elif isinstance(block, str):
                    text_parts.append(block)
            reply_content = "".join(text_parts)

        # --- NEW: Extract and print Token & Tool Stats ---
        in_tokens = final_state.get("total_input_tokens", 0)
        out_tokens = final_state.get("total_output_tokens", 0)
        t_calls = final_state.get("total_tool_calls", 0)
        t_success = final_state.get("successful_tool_calls", 0)

        # Calculate success rate safely to avoid division by zero
        success_rate = (t_success / t_calls * 100) if t_calls > 0 else 0.0

        # Print the formatted metrics summary to the console
        logger.success("AGENT METRICS SUMMARY")
        logger.success(
            f"\n\nInput Tokens: {in_tokens}\nOutput Tokens: {out_tokens}\nTool Calls: {t_success} successful / {t_calls} total\nSuccess Rate:  {success_rate:.1f}%")

        # 4. Return response to the frontend
        return AgentResponse(
            reply=reply_content,
            final_snapshot_id=final_state.get("snapshot_id")
        )

    except Exception as e:
        # Print error to console (Safe now due to UTF-8 fix)
        print(f"[Main Error] {e}")

        # Return a 500 error to the client, but keep the server running
        raise HTTPException(
            status_code=500,
            detail=f"Gardener encountered a technical issue: {str(e)}"
        )


@app.post("/api/config/set_key")
async def set_api_key(request: SetKeyRequest):
    """
    Sets a specific API key securely.
    If is_current is True, updates the default LLM in config and reloads the LLM client.
    """
    try:
        # Validate key type against Enum
        if request.key_type not in APIKey.__members__:
            raise HTTPException(status_code=400, detail=f"Invalid key type: {request.key_type}")

        # Save the secret
        credential_manager.save_secret(APIKey[request.key_type].value, request.key_value)

        # Update current LLM if requested
        if request.is_current:
            with shelve.open(CONFIG_DB) as db:
                db[CUR_LLM] = request.key_type

            # Reload the LLM client to apply changes immediately
            reload_llm()
            logger.info(f"LLM reloaded with new key type: {request.key_type}")

        return {"status": "success", "message": f"Key {request.key_type} saved successfully."}
    except Exception as e:
        logger.error(f"Error setting key: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@app.get("/api/config/check_key")
async def check_api_key(key_type: Optional[str] = None):
    """
    Checks if a specific API key is set.
    If key_type is None, checks the currently configured default LLM.
    Returns boolean status.
    """
    try:
        # If no key provided, try to get from config
        if not key_type:
            with shelve.open(CONFIG_DB) as db:
                key_type = db.get(CUR_LLM)

            if not key_type:
                return {"key_type": None, "is_set": False, "message": "No default LLM configured."}

        if key_type not in APIKey.__members__:
            raise HTTPException(status_code=400, detail=f"Invalid key type: {key_type}")

        secret = credential_manager.get_secret(APIKey[key_type].value)
        is_set = secret is not None and len(secret) > 0
        return {"key_type": key_type, "is_set": is_set}
    except Exception as e:
        logger.error(f"Error checking key: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@app.post("/api/config/verify_key")
async def verify_api_key(request: VerifyKeyRequest):
    """
    Verifies if the API key is valid by making an online check.
    If key_value is provided, it validates that specific key.
    Otherwise, it retrieves the stored key for the given (or default) key_type.
    """
    try:
        key_type = request.key_type
        key_value = request.key_value

        # Case 1: Validate a specific key value provided in the request
        if key_value:
            if not key_type:
                raise HTTPException(status_code=400,
                                    detail="Key type is required when validating a specific key value.")

            if key_type not in APIKey.__members__:
                raise HTTPException(status_code=400, detail=f"Invalid key type: {key_type}")

            # Validate the provided key directly
            is_valid = await validate_api_key_online(APIKey[key_type].value, key_value)
            return {"valid": is_valid, "message": "Key is valid." if is_valid else "Key validation failed."}

        # Case 2: Validate the stored key
        # If no key_type is provided, get the default from config
        if not key_type:
            with shelve.open(CONFIG_DB) as db:
                key_type = db.get(CUR_LLM)
            if not key_type:
                raise HTTPException(status_code=404, detail="No key type provided and no default LLM configured.")

        if key_type not in APIKey.__members__:
            raise HTTPException(status_code=400, detail=f"Invalid key type: {key_type}")

        secret = credential_manager.get_secret(APIKey[key_type].value)
        if not secret:
            return {"valid": False,
                    "message": "The API key is not available in your system’s encrypted credential store."}

        # Perform online validation using stored secret
        is_valid = await validate_api_key_online(APIKey[key_type].value, secret)

        if is_valid:
            return {"valid": True, "message": "Local stored API key is valid."}
        else:
            return {"valid": False, "message": "Stored API key is invalid or failed validation."}

    except HTTPException as he:
        raise he
    except Exception as e:
        logger.error(f"Error verifying key: {e}")
        raise HTTPException(status_code=500, detail=str(e))


if __name__ == "__main__":
    # Start the Uvicorn server on port 41889
    # reload=True is good for dev, consider removing for production
    uvicorn.run("main:app", host="0.0.0.0", port=41889, reload=True)
