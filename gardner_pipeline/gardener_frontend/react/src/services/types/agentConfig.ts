export type AgentApiKeyType =
  | 'OPENAI_API_KEY'
  | 'ANTHROPIC_API_KEY'
  | 'GOOGLE_API_KEY'
  | 'QWEN_API_KEY'
  | 'DEEPSEEK_API_KEY'
  | 'HUGGINGFACE_TOKEN'
  | 'MISTRAL_API_KEY'

export interface SetApiKeyRequest {
  key_type: AgentApiKeyType
  key_value: string
  is_current?: boolean
}

export interface SetApiKeyResponse {
  status: string
  message: string
}

export interface CheckApiKeyResponse {
  key_type: AgentApiKeyType | null
  is_set: boolean
  message?: string
}

export interface VerifyApiKeyRequest {
  key_type?: AgentApiKeyType
  key_value?: string
}

export interface VerifyApiKeyResponse {
  valid: boolean
  message: string
}
