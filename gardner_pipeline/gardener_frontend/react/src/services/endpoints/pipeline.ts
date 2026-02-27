import { API_BASE } from '../config'
import type { RunPipelineRequest, RunPipelineResponse } from '../types'

export async function runPipeline(
  data: RunPipelineRequest
): Promise<RunPipelineResponse> {
  const response = await fetch(`${API_BASE}/pipeline/run`, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(data),
  })

  if (!response.ok) {
    let message = 'Pipeline execution failed'
    try {
      const error = await response.json()
      message = error.message || message
    } catch {
      // keep fallback message
    }
    throw new Error(message)
  }

  return response.json()
}
