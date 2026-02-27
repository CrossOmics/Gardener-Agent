from pathlib import Path
from typing import Optional
from infrastructure.filesystem.constants.filesystem_constants import USER_PROJECT_ROOT


def get_project_relative(project_id: str, subspace: str, filename: Optional[str] = None) -> str:
    """
    Build a POSIX-style relative_key for saving files inside a project workspace.

    Example:
        get_project_relative("projects", "proj123", "qc", "plot.png")
        -> "projects/proj123/qc/plot.png"
    """
    parts = [USER_PROJECT_ROOT, project_id, subspace]
    if filename:
        parts.append(filename)

    return Path(*parts).as_posix()


def get_dataset_relative(
        project_id: str,
        dataset_id: str,
        subspace: str,
        filename: Optional[str] = None
) -> str:
    """
    Build a POSIX-style relative_key for saving files inside a specific dataset folder.

    Structure: USER_PROJECT_ROOT / project_id / dataset_id / subspace / [filename]

    Example:
        get_dataset_relative("p_123", "ds_456", "raw_data", "matrix.h5ad")
        -> "projects/p_123/ds_456/raw_data/matrix.h5ad"
    """
    parts = [USER_PROJECT_ROOT, project_id, dataset_id, subspace]
    if filename:
        parts.append(filename)

    return Path(*parts).as_posix()
