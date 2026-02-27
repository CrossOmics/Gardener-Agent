from typing import Optional, Dict, Any, List
from datetime import datetime
from pydantic import BaseModel, field_validator


class SnapshotResponse(BaseModel):
    id: int
    dataset_id: str

    snapshot_id: str
    branch_name: str
    snapshot_path: str
    parent_snapshot_id: Optional[str] = None
    params_json: Optional[Dict[str, Any]] = None
    thumbnail_json: Optional[Dict[str, Any]] = None
    user_notes: Optional[str] = None
    create_time: datetime
    end_time: Optional[datetime] = None



    class Config:
        from_attributes = True  # Allows Pydantic to read from Peewee models

    @field_validator('dataset_id', mode='before')
    def extract_id_from_object(cls, v: Any) -> str:
        # if v  Dataset Object (Peewee Model)，Take its string ID
        if hasattr(v, 'dataset_id'):
            return str(v.dataset_id)

        return str(v)
