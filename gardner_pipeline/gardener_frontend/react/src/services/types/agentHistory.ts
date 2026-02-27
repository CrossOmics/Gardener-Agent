export interface CreateChatSessionRequest {
  project_id: string
  initial_name?: string | null
  agent_name?: string | null
}

export interface ChatSessionSummary {
  session_id: string
  name: string
  agent?: string
  updated_at?: string
}

export interface CreateChatSessionResponse {
  session_id: string
  name?: string
  status?: string
  created_at?: string
  project_id?: unknown
}

export interface RenameSessionRequest {
  new_name: string
}

export interface CreateAgentMessageRequest {
  message_type: 'user_input' | 'agent_thought' | 'tool_call' | 'tool_result' | 'agent_final'
  data: Record<string, unknown>
}

export interface CreateAgentMessageResponse {
  msg_id: string
  status: string
  type: string
}

export interface SessionHistoryMessage {
  id: string
  type: string
  content: Record<string, unknown>
  timestamp?: string
}

export interface SessionHistoryResponse {
  session_id: string
  total_messages: number
  page: number
  size: number
  messages: SessionHistoryMessage[]
}
