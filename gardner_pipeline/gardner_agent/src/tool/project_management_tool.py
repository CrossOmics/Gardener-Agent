from base.base_call import call_backend
from langchain_core.tools import tool
from typing import Optional


@tool(
    description="Create a new research project workspace. Use this to initialize a new study area.",
    parse_docstring=True
)
async def create_new_project(
        project_name: str,
        description: str = ""
) -> str:
    """
    Initialize a new empty project and a default chat session.

    # CORE RESPONSIBILITY
    Use this tool when the user wants to start a new analysis from scratch or
    organize a new set of datasets into a distinct workspace.

    # CONSTRAINTS
    - Only call this when the user explicitly requests a "new project".
    - Avoid creating duplicate projects with the same name.

    Args:
        project_name: The human-readable name for the new project. REQUIRED.
        description: A brief summary of the project goals or scope.
    """

    try:
        data = await call_backend("POST", "/project/create", {
            "project_name": project_name,
            "description": description
        })

        return (
            f"Project created successfully.\n"
            f"- Project ID: {data.get('project_id')}\n"
            f"- Project Name: {data.get('project_name')}\n"
            f"- Initial Session ID: {data.get('session_id')}\n"
            f"Observation: A new workspace is ready. You can now import datasets to this project."
        )
    except Exception as e:
        return f"Project creation failed: {str(e)}"


@tool(
    description="Get the absolute file system path of the workspace root.",
    parse_docstring=True
)
async def get_workspace_root() -> str:
    """
    Retrieve the base directory where all project data is stored locally.

    # CORE RESPONSIBILITY
    Use this when the user asks where their files are located on the disk or
    when you need to explain the data storage structure.

    # STRATEGY
    This is an informational tool. It does not modify data.
    """
    try:
        root_path = await call_backend("GET", "/project/root")
        return f"The current workspace root is located at: {root_path}"
    except Exception as e:
        return f"Failed to retrieve workspace root: {str(e)}"


@tool(
    description="Permanently delete a dataset and its associated analysis results. Optionally keeps the final snapshot.",
    parse_docstring=True
)
async def delete_dataset(
        dataset_id: str,
        keep_final: bool = False
) -> str:
    """
    Remove a dataset record and all physical files (raw data and snapshots).

    # CORE RESPONSIBILITY
    Use this tool when the researcher wants to discard a specific dataset
    to clean up the project or correct a wrong import.

    # CONSTRAINTS
    - IRREVERSIBLE: Warn the user before calling this if they are unsure.
    - Snapshot Loss: By default, all snapshots derived from this dataset will also be destroyed.

    Args:
        dataset_id: The unique business ID of the dataset to be deleted. REQUIRED.
        keep_final: If true, the latest snapshot of this dataset will be preserved. Defaults to False.
    """
    try:
        # Construct the payload for query parameters
        params = {
            "dataset_id": dataset_id,
            "keep_final": keep_final
        }
        # The endpoint should not contain the dataset_id in the path
        await call_backend("DELETE", "/project/dataset/delete", payload=params)

        if keep_final:
            return f"Dataset {dataset_id} has been successfully removed, but its final snapshot was preserved."
        else:
            return f"Dataset {dataset_id} and all related files have been successfully and permanently removed."

    except Exception as e:
        return f"Dataset deletion failed: {str(e)}"
