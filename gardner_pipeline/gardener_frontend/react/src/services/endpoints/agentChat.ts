import { AGENT_API_BASE } from '../config'
import type { AgentChatRequest, AgentChatResponse } from '../types'

export async function agentChat(data: AgentChatRequest): Promise<AgentChatResponse> {
  const response = await fetch(`${AGENT_API_BASE}/agent/chat`, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(data),
  })

  if (!response.ok) {
    let message = 'Failed to chat with agent'
    try {
      const error = await response.json()
      message = error?.detail || error?.message || message
    } catch {
      // keep fallback
    }
    throw new Error(message)
  }

  return response.json()
}
