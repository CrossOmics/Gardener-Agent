standard_settings = {
    "preprocessing_params": {
        "skip_qc_calculation": False,
        "skip_qc_filter": False,
        "organism": "Human",
        "min_genes": 300,
        "min_cells": 5,
        "pct_mt_max": 15.0,
        "max_counts": 40000,

        "skip_hvg": False,
        "n_top_genes": 2500,
        "flavor": "seurat",
        "target_sum": 10000.0,

        "skip_pca": False,
        "n_comps": 50,
        "svd_solver": "arpack",

        "skip_neighbors": False,
        "n_neighbors": 20,
        "n_pcs": 30
    },

    "clustering_params": {
        "method": "leiden",
        "resolution": 0.5,
        "run_hierarchical": True
    },

    "deg_params": {
        "method": "wilcoxon",
        "n_top_genes": 10
    },

    "annotation_params": {
        "categories": [
            "CellMarker_2024",
            "BioPlanet_2019",
            "NCI-Nature_2016",
            "GO_Biological_Process_2023"
        ],
        "top_n_genes": 50,
        "model_names": ["Immune_All_Low.pkl"],
        "majority_voting": True,
        "target_cluster_col": "leiden",
        "llm_annotation": False,
        "uncertain_threshold": 0.5
    }
}
