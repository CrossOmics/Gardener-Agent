import { API_BASE } from '../config'
import type { PreprocessingFullRequest, PreprocessingFullResponse } from '../types'

export async function runPreprocessingFull(
  data: PreprocessingFullRequest
): Promise<PreprocessingFullResponse> {
  const response = await fetch(`${API_BASE}/preprocessing/run/full`, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(data),
  })

  if (!response.ok) {
    let message = 'Preprocessing failed'
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
