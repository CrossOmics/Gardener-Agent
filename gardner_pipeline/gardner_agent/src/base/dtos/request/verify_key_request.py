from typing import Optional

from pydantic import BaseModel


class VerifyKeyRequest(BaseModel):
    key_type: Optional[str] = None
    key_value: Optional[str] = None

