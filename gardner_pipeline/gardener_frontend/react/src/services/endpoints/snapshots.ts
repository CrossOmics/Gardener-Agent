import { API_BASE } from '../config'
import type {
  DeleteSnapshotsByStageRequest,
  DeleteSnapshotsByStageResponse,
  SnapshotAncestorNode,
  SnapshotBranchName,
  SnapshotQueryResponse,
} from '../types'

const latestStageSnapshotInflight = new Map<string, Promise<SnapshotQueryResponse<Record<string, unknown>, Record<string, unknown>>>>()

async function readErrorMessage(response: Response, fallback: string): Promise<string> {
  try {
    const error = await response.json()
    return error.detail || error.message || fallback
  } catch {
    return fallback
  }
}

export async function getSnapshotDetails<TInput = Record<string, unknown>, TOutput = Record<string, unknown>>(
  snapshotId: string
): Promise<SnapshotQueryResponse<TInput, TOutput>> {
  const response = await fetch(`${API_BASE}/snapshots/query/${encodeURIComponent(snapshotId)}`, {
    method: 'GET',
  })

  if (!response.ok) {
    throw new Error(await readErrorMessage(response, 'Failed to query snapshot details'))
  }

  return response.json()
}

export async function getLatestStageSnapshot<TInput = Record<string, unknown>, TOutput = Record<string, unknown>>(
  datasetId: string,
  branchName: SnapshotBranchName
): Promise<SnapshotQueryResponse<TInput, TOutput>> {
  const requestKey = `${datasetId}::${branchName}`
  const existing = latestStageSnapshotInflight.get(requestKey)
  if (existing) {
    return existing as Promise<SnapshotQueryResponse<TInput, TOutput>>
  }

  const params = new URLSearchParams({
    dataset_id: datasetId,
    branch_name: branchName,
  })
  const request = (async () => {
    const response = await fetch(`${API_BASE}/snapshots/stage/query/?${params.toString()}`, {
      method: 'GET',
    })

    if (!response.ok) {
      throw new Error(await readErrorMessage(response, 'Failed to query latest stage snapshot'))
    }

    return response.json() as Promise<SnapshotQueryResponse<TInput, TOutput>>
  })()

  latestStageSnapshotInflight.set(requestKey, request as Promise<SnapshotQueryResponse<Record<string, unknown>, Record<string, unknown>>>)
  try {
    return await request
  } finally {
    latestStageSnapshotInflight.delete(requestKey)
  }
}

export async function deleteSnapshotsByStage(
  data: DeleteSnapshotsByStageRequest
): Promise<DeleteSnapshotsByStageResponse> {
  const params = new URLSearchParams({
    dataset_id: data.dataset_id,
    stage_name: data.stage_name,
    keep_latest: String(Boolean(data.keep_latest)),
  })

  const response = await fetch(`${API_BASE}/snapshots/stage/delete?${params.toString()}`, {
    method: 'DELETE',
  })

  if (!response.ok) {
    throw new Error(await readErrorMessage(response, 'Failed to delete stage snapshots'))
  }

  return response.json()
}

export async function getSnapshotAncestors(
  snapshotId: string
): Promise<SnapshotAncestorNode[]> {
  const response = await fetch(`${API_BASE}/snapshots/lineage/ancestors/${encodeURIComponent(snapshotId)}`, {
    method: 'GET',
  })

  if (!response.ok) {
    throw new Error(await readErrorMessage(response, 'Failed to query snapshot ancestors'))
  }

  return response.json()
}

export async function deleteSnapshot(
  snapshotId: string
): Promise<{ msg: string }> {
  const response = await fetch(`${API_BASE}/snapshots/delete/${encodeURIComponent(snapshotId)}`, {
    method: 'DELETE',
  })

  if (!response.ok) {
    throw new Error(await readErrorMessage(response, 'Failed to delete snapshot'))
  }

  return response.json()
}
