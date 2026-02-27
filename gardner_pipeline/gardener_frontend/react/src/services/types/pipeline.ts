export interface PreprocessingParams {
  skip_qc_filter?: boolean
  skip_hvg?: boolean
  skip_pca?: boolean
  skip_neighbors?: boolean
  min_genes?: number
  min_cells?: number
  pct_mt_max?: number
  max_counts?: number
  pct_hb_max?: number
  n_top_genes?: number
  flavor?: 'seurat' | 'cell_ranger' | 'seurat_v3'
  target_sum?: number
  n_comps?: number
  svd_solver?: string
  n_neighbors?: number
  n_pcs?: number
}

export interface ClusteringParams {
  method?: 'leiden' | 'louvain' | 'cplearn'
  resolution?: number
  run_hierarchical?: boolean
}

export interface DEGParams {
  groupby?: string
  method?: 'wilcoxon' | 't-test' | 'logreg'
  n_top_genes?: number
}

export interface AnnotationParams {
  categories?: string[]
  top_n_genes?: number
  model_names?: string[]
  majority_voting?: boolean
  target_cluster_col?: string
}

export interface RunPipelineRequest {
  project_id: string
  dataset_id: string
  preprocessing_params?: PreprocessingParams
  clustering_params?: ClusteringParams
  deg_params?: DEGParams
  annotation_params?: AnnotationParams
}

export interface RunPipelineResponse {
  pipeline_run_id: string
  final_snapshot_id: string
  status: string
  msg?: string
  steps_completed: string[]
  snapshots: Record<string, string>
}
