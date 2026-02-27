from pydantic import BaseModel, Field


class RenameSessionRequest(BaseModel):
    new_name: str
