from typing import Dict, Optional
from pydantic import BaseModel


class ClusteringResultDTO(BaseModel):
    snapshot_id: str
    snapshot_path: str

    # Visualization paths
    umap_plot_path: str = None
    dendrogram_path: Optional[str] = None

    # Metadata
    clusters_summary: Dict[str, int] = None
    msg: str = None
