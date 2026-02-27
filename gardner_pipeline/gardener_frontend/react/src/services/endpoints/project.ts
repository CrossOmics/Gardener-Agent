import { API_BASE } from '../config'
import type {
  DisplayNamesRequest,
  DisplayNamesResponse,
  ImportDatasetRequest,
  ImportDatasetResponse,
  ProjectEntityMapResponse,
  RenameEntityRequest,
  RenameEntityResponse,
} from '../types'

const DISPLAY_MAP_CACHE_TTL_MS = 5000

interface ProjectDisplayMapCacheEntry {
  data: ProjectEntityMapResponse
  expiresAt: number
}

const projectDisplayMapCache = new Map<string, ProjectDisplayMapCacheEntry>()
const projectDisplayMapInflight = new Map<string, Promise<ProjectEntityMapResponse>>()

function invalidateProjectDisplayMapCache() {
  projectDisplayMapCache.clear()
  projectDisplayMapInflight.clear()
}

async function readErrorMessage(response: Response, fallback: string): Promise<string> {
  try {
    const error = await response.json()
    return error.detail || error.message || fallback
  } catch {
    return fallback
  }
}

export async function importDataset(
  data: ImportDatasetRequest
): Promise<ImportDatasetResponse> {
  const response = await fetch(`${API_BASE}/project/dataset/import`, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(data),
  })

  if (!response.ok) {
    throw new Error(await readErrorMessage(response, 'Failed to import dataset'))
  }

  return response.json()
}

export async function renameEntityDisplayName(
  data: RenameEntityRequest
): Promise<RenameEntityResponse> {
  const response = await fetch(`${API_BASE}/project/rename`, {
    method: 'PATCH',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(data),
  })

  if (!response.ok) {
    throw new Error(await readErrorMessage(response, 'Failed to rename entity'))
  }

  const payload = await response.json()
  invalidateProjectDisplayMapCache()
  return payload
}

export async function getDisplayNames(
  data: DisplayNamesRequest
): Promise<DisplayNamesResponse> {
  const response = await fetch(`${API_BASE}/project/display/names`, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(data),
  })

  if (!response.ok) {
    throw new Error(await readErrorMessage(response, 'Failed to fetch display names'))
  }

  return response.json()
}

export async function getProjectDisplayMap(
  projectId: string,
  options?: { forceRefresh?: boolean }
): Promise<ProjectEntityMapResponse> {
  const forceRefresh = Boolean(options?.forceRefresh)
  const now = Date.now()
  if (!forceRefresh) {
    const cached = projectDisplayMapCache.get(projectId)
    if (cached && cached.expiresAt > now) {
      return cached.data
    }
    const inflight = projectDisplayMapInflight.get(projectId)
    if (inflight) {
      return inflight
    }
  }

  const request = (async () => {
  const response = await fetch(`${API_BASE}/project/display/map/${encodeURIComponent(projectId)}`, {
    method: 'GET',
  })

  if (!response.ok) {
    throw new Error(await readErrorMessage(response, 'Failed to fetch project display map'))
  }

    const payload = await response.json()
    projectDisplayMapCache.set(projectId, {
      data: payload,
      expiresAt: Date.now() + DISPLAY_MAP_CACHE_TTL_MS,
    })
    return payload
  })()

  projectDisplayMapInflight.set(projectId, request)
  try {
    return await request
  } finally {
    projectDisplayMapInflight.delete(projectId)
  }
}

export async function deleteDataset(
  datasetId: string,
  keepFinal: boolean = false
): Promise<{ msg: string }> {
  const params = new URLSearchParams({
    dataset_id: datasetId,
    keep_final: String(keepFinal),
  })
  const response = await fetch(`${API_BASE}/project/dataset/delete?${params.toString()}`, {
    method: 'DELETE',
  })

  if (!response.ok) {
    throw new Error(await readErrorMessage(response, 'Failed to delete dataset'))
  }

  const payload = await response.json()
  invalidateProjectDisplayMapCache()
  return payload
}
