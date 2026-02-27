from typing import Dict, Any, Optional
from langchain_core.tools import tool
from base.base_call import call_backend


@tool(
    description="Create and save a new user preference profile with specific analysis settings.",
    parse_docstring=True
)
async def create_user_preference(
        settings: Dict[str, Any],
        preference_name: Optional[str] = None
) -> str:
    """
    Save a set of analysis parameters (settings) as a named preference for future use.

    # CORE RESPONSIBILITY
    Use this tool when the user wants to "save" or "remember" a specific configuration
    (e.g., specific QC thresholds or clustering resolutions) so they can reuse it later.

    # STRATEGY
    1. Collect Parameters: Ensure you have the dictionary of settings to save (e.g., {'min_genes': 500, 'resolution': 0.8}).
    2. Naming: If the user provides a name, use it. If not, the system will auto-generate one.

    Args:
        settings: The dictionary of parameters to save (e.g., {"min_genes": 200, "method": "leiden"}). REQUIRED.
        preference_name: A custom name for this preference profile (e.g., "Strict_QC_Config"). Optional.
    """
    try:
        # Construct payload matching CreateUserPreferenceRequest
        payload = {
            "settings": settings,
            "preference_name": preference_name
        }

        data = await call_backend("POST", "/user/preference/create", payload)

        return (
            f"Preference saved successfully.\n"
            f"- Name: {data.get('preference_name')}\n"
            f"- Created At: {data.get('created_at')}\n"
            f"Observation: The settings have been stored. You can load them later using the name '{data.get('preference_name')}'."
        )
    except Exception as e:
        return f"Failed to create user preference: {str(e)}"


@tool(
    description="Retrieve a specific user preference profile by its unique name.",
    parse_docstring=True
)
async def get_user_preference(
        name: str
) -> str:
    """
    Load a saved preference profile to apply its settings.

    # CORE RESPONSIBILITY
    Use this tool when the user asks to "load" or "use" a specific named configuration
    (e.g., "Use the 'Standard_Pipeline' settings").

    # CONSTRAINTS
    - You must know the exact name of the preference.

    Args:
        name: The unique name of the preference to retrieve. REQUIRED.
    """
    try:
        # Pass parameters via query string for GET request
        data = await call_backend("GET", "/user/preference/query", {"name": name})

        settings = data.get("settings", {})
        return (
            f"Preference '{name}' loaded.\n"
            f"- Settings: {settings}\n"
            f"- Last Updated: {data.get('updated_at')}\n"
            f"Observation: You now have the requested parameters. Use these values for the next analysis steps."
        )
    except Exception as e:
        return f"Failed to retrieve preference '{name}': {str(e)}"


@tool(
    description="Fetch the most recently used or created user preference profile.",
    parse_docstring=True
)
async def get_latest_preference() -> str:
    """
    Retrieve the single most recent preference setting.

    # CORE RESPONSIBILITY
    Use this tool when the user says "use my last settings" or "repeat with previous config"
    without specifying a name.

    # STRATEGY
    - This automatically fetches the record with the latest timestamp.
    """
    try:
        data = await call_backend("GET", "/user/preference/latest/one")

        name = data.get("preference_name")
        settings = data.get("settings", {})

        return (
            f"Latest preference found: '{name}'.\n"
            f"- Settings: {settings}\n"
            f"Observation: These are the most recently used settings. You can apply them now."
        )
    except Exception as e:
        return f"Failed to fetch latest preference: {str(e)}"


@tool(
    description="Update an existing user preference with new settings or rename it.",
    parse_docstring=True
)
async def update_user_preference(
        name: str,
        settings: Dict[str, Any],
        new_name: Optional[str] = None
) -> str:
    """
    Modify an existing preference profile.

    # CORE RESPONSIBILITY
    Use this tool when the user wants to change the parameters stored in a saved profile
    or rename it.

    Args:
        name: The CURRENT unique name of the preference to update. REQUIRED.
        settings: The NEW dictionary of parameters to overwrite the old ones. REQUIRED.
        new_name: A new name for the preference (if renaming is requested). Optional.
    """
    try:
        # Construct payload matching UpdateUserPreferenceRequest
        payload = {
            "name": name,
            "settings": settings,
            "new_name": new_name
        }

        data = await call_backend("PATCH", "/user/preference/update", payload)

        return (
            f"Preference '{name}' updated successfully.\n"
            f"- Current Name: {data.get('preference_name')}\n"
            f"- New Settings: {data.get('settings')}\n"
            f"Observation: The profile has been updated."
        )
    except Exception as e:
        return f"Failed to update preference '{name}': {str(e)}"
