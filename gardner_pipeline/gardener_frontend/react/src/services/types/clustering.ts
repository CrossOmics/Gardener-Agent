export interface ClusteringCreateRequest {
  project_id: string
  dataset_id: string
  snapshot_id?: string
  method?: 'leiden' | 'louvain' | 'cplearn'
  resolution?: number
  run_hierarchical?: boolean
}

export interface ClusteringCreateResponse {
  snapshot_id: string
  snapshot_path: string
  umap_plot_path?: string
  dendrogram_path?: string | null
  clusters_summary?: Record<string, number>
  msg: string
}

export interface ClusteringMergeRequest {
  project_id: string
  dataset_id: string
  snapshot_id: string
  method?: 'leiden' | 'louvain' | 'cplearn'
  clusters_to_merge: Array<string | number>
  new_label?: string
  overwrite?: boolean
}

export interface ClusteringMergeResponse {
  snapshot_id: string
  snapshot_path: string
  msg: string
}

export interface ClusteringSubclusterRequest {
  project_id: string
  dataset_id: string
  snapshot_id: string
  source_cluster_col: string
  target_clusters: Array<string | number>
  method?: 'leiden' | 'louvain' | 'cplearn'
  resolution?: number
  overwrite?: boolean
}

export interface ClusteringSubclusterResponse {
  snapshot_id: string
  snapshot_path: string
  msg: string
}
