export interface ImportDatasetRequest {
  project_id: string
  local_file_path: string
  dataset_name?: string
  organism?: string
  tissue_type?: string
  description?: string
}

export interface ImportDatasetResponse {
  project_id: string
  project_name: string
  dataset_id: string
  dataset_name: string
  workspace_path: string
}

export interface CreateProjectRequest {
  project_name: string
  description?: string
}

export interface CreateProjectResponse {
  project_id: string
  project_name: string
  dataset_id: string
  dataset_name: string
  workspace_path: string
}

export type EntityIdType = 'project' | 'dataset' | 'snapshot'

export interface RenameEntityRequest {
  id_type: EntityIdType
  current_id: string
  new_name: string
}

export interface RenameEntityResponse {
  msg: string
  id_type: EntityIdType
  current_id: string
  new_name: string
}

export interface DisplayNamesRequest {
  id_type: EntityIdType
  current_ids: string[]
}

export interface DisplayNamesResponse {
  id_type: EntityIdType
  names: Record<string, string>
}

export interface DatasetEntityMap {
  dataset_name: string
  snapshots: Record<string, Record<string, string>>
}

export interface ProjectEntityMapResponse {
  project_id: string
  project_name: string
  datasets: Record<string, DatasetEntityMap>
}
