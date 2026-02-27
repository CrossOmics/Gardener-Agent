export interface AgentChatRequest {
  message: string
  snapshot_id: string
  project_id: string
  dataset_id?: string | null
  session_id?: string | null
}

export interface AgentChatResponse {
  reply: string
  final_snapshot_id?: string | null
}
