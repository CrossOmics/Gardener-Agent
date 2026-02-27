export interface DGERequest {
  project_id: string
  dataset_id: string
  snapshot_id?: string
  groupby?: string
  method?: 'wilcoxon' | 't-test' | 'logreg'
  n_top_genes?: number
  use_raw?: boolean
}

export interface DGEResponse {
  snapshot_id: string
  snapshot_path: string
  csv_path?: string
  results_file?: string
  rank_genes_plot?: string
  dotplot?: string
  heatmap?: string
  violin?: string
  top_markers?: Record<string, string[]>
  msg: string
}
