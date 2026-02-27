export interface AnnotationFullRequest {
  project_id: string
  dataset_id: string
  snapshot_id: string
  categories?: string[]
  top_n_genes?: number
  model_names?: string[]
  majority_voting?: boolean
  target_cluster_col?: string
}

export interface AnnotationFullResponse {
  snapshot_id: string
  snapshot_path: string
  cluster_id: string | string[] | null
  msg: string
  enrichment_results?: {
    report_path?: string
    result_keys?: string[]
    plots?: string[]
  } & Record<string, unknown>
}

export interface AnnotationSearchRequest {
  keyword?: string
  type?: 'all' | 'celltypist' | 'gseapy'
}

export interface AnnotationSearchItem {
  id: number | string
  method_name: string
  type: 'celltypist' | 'gseapy'
  description?: string
}

export type AnnotationSearchResponse = AnnotationSearchItem[]

export interface AnnotationSelectionSaveRequest {
  project_id: string
  dataset_id: string
  selected_model_ids: string[]
}

export interface AnnotationSelectionSaveResponse {
  msg: string
}

export interface AnnotationLabelsUpdateRequest {
  snapshot_id: string
  file_name: string
  updated_annotation: Record<string, string>
}

export interface AnnotationLabelsUpdateResponse {
  msg: string
  snapshot_id: string
  updated_count: number
  target_file: string
}
