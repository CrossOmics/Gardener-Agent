import json

import httpx
from typing import Any, Dict, Optional
from loguru import logger

from config import BACKEND_BASE_URL


async def call_backend(method: str, endpoint: str, payload: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    """
    Generic async HTTP client for backend communication.

    CRITICAL FEATURES:
    1. Smart Parameter Handling: Automatically puts payload into 'params' (Query) or 'json' (Body) based on method.
    2. Rich Error Extraction: Parses FastAPI 'detail' so the Agent understands the specific failure reason.
    """
    url = f"{BACKEND_BASE_URL}{endpoint}"
    method_upper = method.upper()

    # Smart Assignment
    # GET/DELETE: payload -> URL parameters (e.g., ?id=123)
    # POST/PATCH/PUT: payload -> JSON body
    params = None
    json_body = None

    if method_upper in ["GET", "DELETE"]:
        params = payload
    else:
        json_body = payload

    async with httpx.AsyncClient(timeout=100000.0) as client:
        try:
            # Execute request with correctly mapped data
            response = await client.request(method, url, json=json_body, params=params)

            # --- 2. Rich Error Handling ---
            if response.is_error:
                error_header = f"HTTP {response.status_code}"

                try:
                    # Attempt to parse FastAPI JSON error response
                    error_data = response.json()

                    if "detail" in error_data:
                        detailed_msg = error_data["detail"]

                        # Handle Pydantic validation errors (which return a list)
                        if isinstance(detailed_msg, list):
                            detailed_msg = json.dumps(detailed_msg)

                        error_msg = f"{error_header}: {detailed_msg}"
                    else:
                        # Fallback: JSON exists but no 'detail' key
                        error_msg = f"{error_header}: {response.text}"

                except json.JSONDecodeError:
                    # Fallback: Response is not JSON (e.g., raw 500 HTML or text)
                    error_msg = f"{error_header}: {response.text}"

                logger.error(f"[Backend Error] {error_msg}")

                # Raise exception with the specific message for the Agent's "Reflection" phase
                raise Exception(error_msg)

            # Success
            return response.json()

        except httpx.RequestError as e:
            # Handle low-level network issues (DNS, Connection Refused)
            logger.error(f"Network Error connecting to {url}: {str(e)}")
            raise Exception(f"Network Connection Failed: {str(e)}")

        except Exception as e:
            # Re-raise the detailed exception created above
            raise e
