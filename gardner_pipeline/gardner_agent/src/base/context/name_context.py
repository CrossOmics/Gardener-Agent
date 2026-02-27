from typing import List, Dict
from base.base_call import call_backend
from loguru import logger


async def get_display_names(ids: List[str], id_type: str) -> Dict[str, str]:
    """
    Retrieves the user-friendly display names for a list of entity IDs.

    Args:
        ids (List[str]): A list of unique identifiers (e.g., snapshot IDs, dataset IDs).
        id_type (str): The type of entity. Must be one of 'project', 'dataset', or 'snapshot'.

    Returns:
        Dict[str, str]: A dictionary mapping each ID to its corresponding display name.
                        If a name cannot be found, the value will be an empty string.
    """
    if not ids:
        return {}

    payload = {
        "id_type": id_type,
        "current_ids": ids
    }

    try:
        # Call the backend endpoint defined in the controller
        response = await call_backend("POST", "/project/display/names", payload=payload)
        return response.get("names", {})
    except Exception as e:
        logger.error(f"Failed to fetch display names for {id_type}: {e}")
        return {}
