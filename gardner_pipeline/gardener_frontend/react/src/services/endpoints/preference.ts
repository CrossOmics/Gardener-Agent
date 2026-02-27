import { API_BASE } from '../config'
import type { UpdateUserPreferenceRequest, UserPreferenceResponse } from '../types'

export async function getLatestUserPreference(): Promise<UserPreferenceResponse> {
  const response = await fetch(`${API_BASE}/user/preference/latest/one`, {
    method: 'GET',
  })

  if (!response.ok) {
    let message = 'Failed to load latest user preference'
    try {
      const error = await response.json()
      message = error.detail || error.message || message
    } catch {
      // keep fallback message
    }
    throw new Error(message)
  }

  return response.json()
}

export async function updateUserPreference(
  data: UpdateUserPreferenceRequest
): Promise<UserPreferenceResponse> {
  const response = await fetch(`${API_BASE}/user/preference/update`, {
    method: 'PATCH',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(data),
  })

  if (!response.ok) {
    let message = 'Failed to update user preference'
    try {
      const error = await response.json()
      message = error.detail || error.message || message
    } catch {
      // keep fallback message
    }
    throw new Error(message)
  }

  return response.json()
}
