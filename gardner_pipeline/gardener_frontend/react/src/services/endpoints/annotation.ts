import { API_BASE } from '../config'
import type {
  AnnotationFullRequest,
  AnnotationFullResponse,
  AnnotationLabelsUpdateRequest,
  AnnotationLabelsUpdateResponse,
  AnnotationSearchRequest,
  AnnotationSearchResponse,
  AnnotationSelectionSaveRequest,
  AnnotationSelectionSaveResponse,
} from '../types'

export async function runAnnotationFull(
  data: AnnotationFullRequest
): Promise<AnnotationFullResponse> {
  const response = await fetch(`${API_BASE}/annotation/full`, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(data),
  })

  if (!response.ok) {
    let message = 'Annotation full run failed'
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

export async function searchAnnotationOptions(
  data: AnnotationSearchRequest = { type: 'all' }
): Promise<AnnotationSearchResponse> {
  const params = new URLSearchParams()
  if (data.type) {
    params.set('type', data.type)
  }
  if (data.keyword && data.keyword.trim() !== '') {
    params.set('keyword', data.keyword.trim())
  }
  const query = params.toString()
  const response = await fetch(`${API_BASE}/methods/annotation/search${query ? `?${query}` : ''}`, {
    method: 'GET',
  })

  if (!response.ok) {
    let message = 'Annotation search failed'
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

export async function saveAnnotationSelection(
  data: AnnotationSelectionSaveRequest
): Promise<AnnotationSelectionSaveResponse> {
  const response = await fetch(`${API_BASE}/annotation/selection/save`, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(data),
  })

  if (!response.ok) {
    let message = 'Save annotation selection failed'
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

export async function updateAnnotationLabels(
  data: AnnotationLabelsUpdateRequest
): Promise<AnnotationLabelsUpdateResponse> {
  const response = await fetch(`${API_BASE}/annotation/labels/update`, {
    method: 'PATCH',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(data),
  })

  if (!response.ok) {
    let message = 'Update annotation labels failed'
    try {
      const error = await response.json()
      message = error?.message || error?.detail || message
    } catch {
      // keep fallback
    }
    throw new Error(message)
  }

  return response.json()
}
