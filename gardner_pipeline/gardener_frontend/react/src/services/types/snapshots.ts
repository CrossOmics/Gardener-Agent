export type SnapshotBranchName = 'Preprocessing' | 'Clustering' | 'DGE' | 'Annotation'

export interface SnapshotQueryResponse<TInput = Record<string, unknown>, TOutput = Record<string, unknown>> {
  snapshot_id: string
  parent_snapshot_id?: string | null
  branch_name: SnapshotBranchName
  create_time: string
  user_notes?: string
  input: TInput
  output: TOutput
}

export interface DeleteSnapshotsByStageRequest {
  dataset_id: string
  stage_name: string
  keep_latest?: boolean
}

export interface DeleteSnapshotsByStageResponse {
  msg: string
  deleted_count: number
}

export interface SnapshotAncestorNode {
  snapshot_id: string
  snapshot_name: string
  stage_name: string
  create_time: string
}
