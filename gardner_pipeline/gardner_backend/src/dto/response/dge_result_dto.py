from typing import Dict, List, Optional
from pydantic import BaseModel


class DGEResultDTO(BaseModel):
    snapshot_id: str
    snapshot_path: str

    # Data Export
    csv_path: str

    # Visualization Paths
    rank_genes_plot: str
    dotplot: str
    heatmap: str
    violin: str

    # Lightweight Preview (Top 5 genes per cluster)
    # Example: { "0": ["GeneA", "GeneB"], "1": ["GeneC", "GeneD"] }
    top_markers: Dict[str, List[str]]

    msg: str