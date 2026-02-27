import { API_BASE } from '../config'
import type {
  ChatSessionSummary,
  CreateAgentMessageRequest,
  CreateAgentMessageResponse,
  CreateChatSessionRequest,
  CreateChatSessionResponse,
  RenameSessionRequest,
  SessionHistoryMessage,
  SessionHistoryResponse,
} from '../types'

function normalizeSession(item: unknown): ChatSessionSummary | null {
  if (!item || typeof item !== 'object') return null
  const obj = item as Record<string, unknown>
  const sessionId = String(obj.session_id ?? '').trim()
  if (!sessionId) return null
  return {
    session_id: sessionId,
    name: String(obj.name ?? 'New Bio Analysis').trim() || 'New Bio Analysis',
    agent: typeof obj.agent === 'string' ? obj.agent : undefined,
    updated_at: typeof obj.updated_at === 'string' ? obj.updated_at : undefined,
  }
}

export async function listProjectSessions(projectId: string): Promise<ChatSessionSummary[]> {
  const params = new URLSearchParams({ project_id: projectId })
  const response = await fetch(`${API_BASE}/agent/history/project/session?${params.toString()}`, {
    method: 'GET',
  })

  if (!response.ok) {
    let message = 'Failed to load session list'
    try {
      const error = await response.json()
      message = error.detail || error.message || message
    } catch {
      // keep fallback
    }
    throw new Error(message)
  }

  const payload = await response.json() as unknown
  const rawList = Array.isArray(payload)
    ? payload
    : (payload && typeof payload === 'object' && Array.isArray((payload as { value?: unknown[] }).value))
      ? (payload as { value: unknown[] }).value
      : []

  return rawList
    .map(normalizeSession)
    .filter((x): x is ChatSessionSummary => Boolean(x))
}

export async function createChatSession(
  data: CreateChatSessionRequest
): Promise<CreateChatSessionResponse> {
  const response = await fetch(`${API_BASE}/agent/history/sessions/create`, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(data),
  })

  if (!response.ok) {
    let message = 'Failed to create chat session'
    try {
      const error = await response.json()
      message = error.detail || error.message || message
    } catch {
      // keep fallback
    }
    throw new Error(message)
  }

  return response.json()
}

export async function renameChatSession(
  sessionId: string,
  data: RenameSessionRequest
): Promise<void> {
  const response = await fetch(`${API_BASE}/agent/history/sessions/rename/${sessionId}`, {
    method: 'PATCH',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(data),
  })

  if (!response.ok) {
    let message = 'Failed to rename chat session'
    try {
      const error = await response.json()
      message = error.detail || error.message || message
    } catch {
      // keep fallback
    }
    throw new Error(message)
  }
}

export async function deleteChatSession(sessionId: string): Promise<void> {
  const response = await fetch(`${API_BASE}/agent/history/sessions/delete/${sessionId}`, {
    method: 'DELETE',
  })

  if (!response.ok) {
    let message = 'Failed to delete chat session'
    try {
      const error = await response.json()
      message = error.detail || error.message || message
    } catch {
      // keep fallback
    }
    throw new Error(message)
  }
}

function normalizeHistoryMessage(item: unknown): SessionHistoryMessage | null {
  if (!item || typeof item !== 'object') return null
  const obj = item as Record<string, unknown>
  const id = String(obj.id ?? '').trim()
  if (!id) return null
  return {
    id,
    type: String(obj.type ?? 'user_input'),
    content: (obj.content && typeof obj.content === 'object')
      ? obj.content as Record<string, unknown>
      : {},
    timestamp: typeof obj.timestamp === 'string' ? obj.timestamp : undefined,
  }
}

export async function getSessionHistory(
  sessionId: string,
  page: number = 1,
  size: number = 50
): Promise<SessionHistoryResponse> {
  const params = new URLSearchParams({
    page: String(page),
    size: String(size),
  })
  const response = await fetch(`${API_BASE}/agent/history/${sessionId}/history/all?${params.toString()}`, {
    method: 'GET',
  })

  if (!response.ok) {
    let message = 'Failed to load session history'
    try {
      const error = await response.json()
      message = error.detail || error.message || message
    } catch {
      // keep fallback
    }
    throw new Error(message)
  }

  const payload = await response.json() as {
    session_id?: string
    total_messages?: number
    page?: number
    size?: number
    messages?: unknown[]
  }

  return {
    session_id: payload.session_id || sessionId,
    total_messages: typeof payload.total_messages === 'number' ? payload.total_messages : 0,
    page: typeof payload.page === 'number' ? payload.page : page,
    size: typeof payload.size === 'number' ? payload.size : size,
    messages: Array.isArray(payload.messages)
      ? payload.messages.map(normalizeHistoryMessage).filter((x): x is SessionHistoryMessage => Boolean(x))
      : [],
  }
}

export async function createSessionMessage(
  sessionId: string,
  data: CreateAgentMessageRequest
): Promise<CreateAgentMessageResponse> {
  const response = await fetch(`${API_BASE}/agent/history/${sessionId}/message/create`, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(data),
  })

  if (!response.ok) {
    let message = 'Failed to create session message'
    try {
      const error = await response.json()
      message = error.detail || error.message || message
    } catch {
      // keep fallback
    }
    throw new Error(message)
  }

  return response.json()
}
