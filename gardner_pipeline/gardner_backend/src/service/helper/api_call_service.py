import json
import httpx
from typing import Any, Dict, Optional
from loguru import logger

# Temporary configuration for Agent Server
AGENT_BASE_URL = "http://localhost:41889/api"


class AgentAPIService:
    """
    Service to communicate with the Gardener Agent Server.
    Mimics the logic of base_call.py in the agent.
    """

    @staticmethod
    async def call_agent(method: str, endpoint: str, payload: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
        """
        Generic async HTTP client for Agent communication.
        
        Args:
            method: HTTP method (GET, POST, etc.)
            endpoint: API endpoint (e.g., '/api/agent/chat')
            payload: Data to send (query params for GET/DELETE, JSON body for others)
            
        Returns:
            Dict containing the JSON response from the agent.
        """
        url = f"{AGENT_BASE_URL}{endpoint}"
        method_upper = method.upper()

        # Smart Assignment
        # GET/DELETE: payload -> URL parameters
        # POST/PATCH/PUT: payload -> JSON body
        params = None
        json_body = None

        if method_upper in ["GET", "DELETE"]:
            params = payload
        else:
            json_body = payload

        # High timeout for LLM operations
        async with httpx.AsyncClient(timeout=300.0) as client:
            try:
                logger.info(f"Calling Agent: {method} {url}")
                response = await client.request(method, url, json=json_body, params=params)

                # Error Handling
                if response.is_error:
                    error_header = f"HTTP {response.status_code}"
                    try:
                        error_data = response.json()
                        if "detail" in error_data:
                            detailed_msg = error_data["detail"]
                            if isinstance(detailed_msg, list):
                                detailed_msg = json.dumps(detailed_msg)
                            error_msg = f"{error_header}: {detailed_msg}"
                        else:
                            error_msg = f"{error_header}: {response.text}"
                    except json.JSONDecodeError:
                        error_msg = f"{error_header}: {response.text}"

                    logger.error(f"[Agent Error] {error_msg}")
                    raise Exception(error_msg)

                return response.json()

            except httpx.RequestError as e:
                logger.error(f"Network Error connecting to {url}: {str(e)}")
                raise Exception(f"Network Connection Failed: {str(e)}")
            except Exception as e:
                logger.error(f"Error calling agent: {str(e)}")
                raise e
