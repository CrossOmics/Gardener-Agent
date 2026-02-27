import path from 'path'

export const MAIN_API_BASE_URL = 'http://127.0.0.1:41888'
export const API_BASE = `${MAIN_API_BASE_URL}/api/v1`

export function getPreloadPath(baseDir: string): string {
  return path.join(baseDir, 'preload.js')
}

function extractDataRoot(payload: unknown): string | null {
  if (typeof payload === 'string' && payload.trim() !== '') {
    return payload.trim()
  }
  if (!payload || typeof payload !== 'object') {
    return null
  }
  const record = payload as Record<string, unknown>
  const candidates = [
    record.data_root,
    record.root_path,
    record.project_root,
    record.path,
    record.value,
  ]
  for (const item of candidates) {
    if (typeof item === 'string' && item.trim() !== '') {
      return item.trim()
    }
  }
  return null
}

export async function resolveDataRoot(): Promise<string> {
  const response = await fetch(`${API_BASE}/project/root`, { method: 'GET' })
  if (!response.ok) {
    throw new Error(`Failed to load project root: HTTP ${response.status}`)
  }
  const payload = await response.json()
  const resolved = extractDataRoot(payload)
  if (!resolved) {
    throw new Error('Failed to load project root: invalid response payload')
  }
  return resolved
}

export function getPaths(dataRoot: string, baseDir: string) {
  return {
    DATA_ROOT: dataRoot,
    PROJECTS_PATH: path.join(dataRoot, 'projects'),
    PRELOAD_PATH: getPreloadPath(baseDir),
  }
}
