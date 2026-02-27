export interface PreprocessingFullRequest {
  project_id: string
  dataset_id: string
  organism?: string
  custom_prefixes?: Record<string, string> | null
  skip_qc_calculation?: boolean
  min_genes?: number
  cell_min_counts?: number
  max_genes?: number
  cell_max_counts?: number
  min_cells?: number
  gene_min_counts?: number
  max_cells?: number
  gene_max_counts?: number
  pct_mt_max?: number
  pct_hb_max?: number
  skip_qc_filter?: boolean
  n_top_genes?: number
  flavor?: 'seurat' | 'cell_ranger' | 'seurat_v3'
  target_sum?: number
  skip_hvg?: boolean
  n_comps?: number
  svd_solver?: string
  skip_pca?: boolean
  n_neighbors?: number
  n_pcs?: number
  skip_neighbors?: boolean
}

export interface PreprocessingFullResponse {
  project_id?: string
  dataset_id?: string
  snapshot_id?: string
  snapshot_path?: string
  final_stage?: string
  executed_steps?: string[]
  skipped_steps?: string[]
  anndata_cache_path?: string
  message?: string
  msg?: string
}
