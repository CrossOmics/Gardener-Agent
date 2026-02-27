from pydantic import BaseModel


class RenameResponse(BaseModel):
    """
    Data Transfer Object for the rename operation response.
    """
    msg: str
    id_type: str
    current_id: str
    new_name: str
