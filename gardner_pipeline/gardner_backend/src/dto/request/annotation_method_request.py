from typing import Optional, List
from pydantic import BaseModel, Field


class SearchAnnotationMethodRequest(BaseModel):
    """
    Request DTO for searching annotation reference methods.
    """
    keyword: Optional[str] = Field(default="", description="Search term for name, description, species, or organ.")
    type: str = Field(default="all", description="Filter by type: 'gseapy', 'celltypist', or 'all'.")
