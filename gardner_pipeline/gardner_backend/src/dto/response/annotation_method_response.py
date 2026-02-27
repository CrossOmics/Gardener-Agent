from typing import Optional
from pydantic import BaseModel
class AnnotationMethodResponse(BaseModel):
    """
    Response DTO returning method details.
    """
    id: int
    method_name: str
    type: str
    description: Optional[str] = None
    species: str
    organ: Optional[str] = None

    class Config:
        # Pydantic V2 configuration to allow creation from ORM objects
        from_attributes = True
