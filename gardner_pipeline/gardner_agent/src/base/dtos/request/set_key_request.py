from pydantic import BaseModel


class SetKeyRequest(BaseModel):
    key_type: str
    key_value: str = ""
    is_current: bool = True