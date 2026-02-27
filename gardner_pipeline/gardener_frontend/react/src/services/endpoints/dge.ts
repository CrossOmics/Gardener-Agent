import { API_BASE } from '../config'
import type { DGERequest, DGEResponse } from '../types'

export async function runDGE(
  data: DGERequest
): Promise<DGEResponse> {
  const response = await fetch(`${API_BASE}/dge/new`, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(data),
  })

  if (!response.ok) {
    let message = 'DGE analysis failed'
    try {
      const error = await response.json()
      message = error.message || message
    } catch {
      // keep fallback
    }
    throw new Error(message)
  }

  return response.json()
}
