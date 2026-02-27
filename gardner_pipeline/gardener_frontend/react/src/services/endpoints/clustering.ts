import { API_BASE } from '../config'
import type {
  ClusteringCreateRequest,
  ClusteringCreateResponse,
  ClusteringMergeRequest,
  ClusteringMergeResponse,
  ClusteringSubclusterRequest,
  ClusteringSubclusterResponse,
} from '../types'

export async function runClusteringCreate(
  data: ClusteringCreateRequest
): Promise<ClusteringCreateResponse> {
  const response = await fetch(`${API_BASE}/clustering/create`, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(data),
  })

  if (!response.ok) {
    let message = 'Clustering creation failed'
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

export async function runClusteringMerge(
  data: ClusteringMergeRequest
): Promise<ClusteringMergeResponse> {
  const response = await fetch(`${API_BASE}/clustering/merge`, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(data),
  })

  if (!response.ok) {
    let message = 'Cluster merge failed'
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

export async function runClusteringSubcluster(
  data: ClusteringSubclusterRequest
): Promise<ClusteringSubclusterResponse> {
  const response = await fetch(`${API_BASE}/clustering/sub`, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(data),
  })

  if (!response.ok) {
    let message = 'Sub-clustering failed'
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
