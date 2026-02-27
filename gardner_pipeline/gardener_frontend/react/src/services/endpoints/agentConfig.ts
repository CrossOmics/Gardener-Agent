import { AGENT_API_BASE } from '../config'
import type {
  CheckApiKeyResponse,
  SetApiKeyRequest,
  SetApiKeyResponse,
  VerifyApiKeyRequest,
  VerifyApiKeyResponse,
} from '../types'

async function parseError(response: Response, fallback: string): Promise<string> {
  try {
    const payload = await response.json()
    return payload.detail || payload.message || fallback
  } catch {
    return fallback
  }
}

export async function setApiKey(data: SetApiKeyRequest): Promise<SetApiKeyResponse> {
  const response = await fetch(`${AGENT_API_BASE}/config/set_key`, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(data),
  })

  if (!response.ok) {
    throw new Error(await parseError(response, 'Failed to set API key'))
  }
  return response.json()
}

export async function checkApiKey(keyType?: string): Promise<CheckApiKeyResponse> {
  const params = new URLSearchParams()
  if (keyType) params.set('key_type', keyType)
  const suffix = params.toString() ? `?${params.toString()}` : ''

  const response = await fetch(`${AGENT_API_BASE}/config/check_key${suffix}`, {
    method: 'GET',
  })

  if (!response.ok) {
    throw new Error(await parseError(response, 'Failed to check API key'))
  }
  return response.json()
}

export async function verifyApiKey(data: VerifyApiKeyRequest): Promise<VerifyApiKeyResponse> {
  const response = await fetch(`${AGENT_API_BASE}/config/verify_key`, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(data),
  })

  if (!response.ok) {
    throw new Error(await parseError(response, 'Failed to verify API key'))
  }
  return response.json()
}

