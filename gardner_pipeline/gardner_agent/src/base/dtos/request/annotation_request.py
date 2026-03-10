from typing import Dict, List, Any, Optional
from pydantic import BaseModel

class AnnotationRequest(BaseModel):
    """
    Request payload for LLM-based annotation refinement.
    """
    uncertain_data: Dict[str, Dict[str, Any]]
    valid_cell_types: Optional[str] = None
    # Structure example:
    # {
    #   "9": {
    #       "dge": ["GeneA", "GeneB", ...],
    #       "uncertain_labels": ["LabelX", "LabelY"]
    #   }
    # }
