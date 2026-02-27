from typing import Optional

from pydantic import BaseModel


class ClusterParameter(BaseModel):
    '''
    Construct the clustering call input parameters
    '''
    method: str
    resolution: float
    # Options
    run_hierarchical: bool
    # cplearn parapmeters
    qnn_size: Optional[int] = None
    neighbor_radius: Optional[int] = None
