export interface UserPreferenceResponse {
  id: number
  preference_name: string
  settings: Record<string, unknown>
  created_at: string
  updated_at: string
}

export interface UpdateUserPreferenceRequest {
  name: string
  settings: Record<string, unknown>
  new_name?: string | null
}
