from typing import Optional, Dict, Any, List
from pydantic import BaseModel


class CreateSnapshotRequest(BaseModel):
    dataset_id: str
    snapshot_id: str
    branch_name: str
    snapshot_path: str
    parent_snapshot_id: Optional[str] = None
    params_json: Optional[Dict[str, Any]] = {}
    thumbnail_json: Optional[Dict[str, Any]] = {}
    user_notes: Optional[str] = None


class UpdateSnapshotRequest(BaseModel):
    branch_name: Optional[str] = None
    user_notes: Optional[str] = None
    params_update: Optional[Dict[str, Any]] = None
