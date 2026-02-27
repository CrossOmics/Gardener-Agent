import { useState, useEffect, useCallback } from 'react'
import { Button, Dropdown, Input, Input as AntInput, Modal, Select, message } from 'antd'
import { SendOutlined, CloseOutlined, EllipsisOutlined } from '@ant-design/icons'
import {
  listProjectSessions,
  createChatSession,
  renameChatSession,
  deleteChatSession,
  getSessionHistory,
  agentChat,
  setApiKey,
  verifyApiKey,
} from '../../services/endpoints'
import type { AgentApiKeyType, ChatSessionSummary } from '../../services/types'
import SessionToolbar from './SessionToolbar'
import ChatMessageContent from './ChatMessageContent'
import './styles.css'

const { TextArea } = Input

interface ChatPanelProps {
  projectId?: string
  datasetId?: string | null
  snapshotId?: string | null
  showHeader?: boolean
  onClose?: () => void
  onAgentResponse?: () => void
}

interface Message {
  id: string
  role: 'user' | 'assistant'
  content: string
  pending?: boolean
}

function normalizeAgentReply(reply: string): string {
  const trimmed = (reply || '').trim()
  if (!trimmed) return ''

  let content = trimmed

  // Handle payloads wrapped as JSON string literals, e.g. "\"## Title\\n...\""
  if (
    (content.startsWith('"') && content.endsWith('"')) ||
    (content.startsWith("'") && content.endsWith("'"))
  ) {
    try {
      content = JSON.parse(content)
    } catch {
      content = content.slice(1, -1)
    }
  }

  if (content.startsWith('{') && content.endsWith('}')) {
    try {
      const parsed = JSON.parse(content) as Record<string, unknown>
      if (typeof parsed.report === 'string' && parsed.report.trim() !== '') {
        content = parsed.report
      } else if (typeof parsed.reply === 'string' && parsed.reply.trim() !== '') {
        content = parsed.reply
      }
    } catch {
      // keep original content
    }
  }

  return content
    .replace(/\\n/g, '\n')
    .replace(/\\"/g, '"')
}

function pickStringField(obj: Record<string, unknown>, keys: string[]): string | null {
  for (const key of keys) {
    const value = obj[key]
    if (typeof value === 'string' && value.trim() !== '') {
      return value
    }
  }
  return null
}

function extractHistoryMessageText(content: Record<string, unknown>): string {
  const direct = pickStringField(content, ['content', 'text', 'report', 'reply', 'message'])
  if (direct) return direct

  const nestedContent = content.content
  if (nestedContent && typeof nestedContent === 'object') {
    const nested = nestedContent as Record<string, unknown>
    const nestedText = pickStringField(nested, ['content', 'text', 'report', 'reply', 'message'])
    if (nestedText) return nestedText
  }

  return ''
}

const MODEL_KEY_OPTIONS: AgentApiKeyType[] = [
  'OPENAI_API_KEY',
  'ANTHROPIC_API_KEY',
  'GOOGLE_API_KEY',
  'QWEN_API_KEY',
  'DEEPSEEK_API_KEY',
  'HUGGINGFACE_TOKEN',
  'MISTRAL_API_KEY',
]

function getModelKeyDisplayLabel(keyType: AgentApiKeyType): string {
  return keyType.replace(/_API_KEY$/, '')
}

export default function ChatPanel({
  projectId,
  datasetId,
  snapshotId,
  showHeader = false,
  onClose,
  onAgentResponse,
}: ChatPanelProps) {
  const [inputValue, setInputValue] = useState('')
  const [messages, setMessages] = useState<Message[]>([])
  const [sessions, setSessions] = useState<ChatSessionSummary[]>([])
  const [activeSessionId, setActiveSessionId] = useState<string>('')
  const [loadingSessions, setLoadingSessions] = useState(false)
  const [creatingSession, setCreatingSession] = useState(false)
  const [renamingSession, setRenamingSession] = useState(false)
  const [deletingSession, setDeletingSession] = useState(false)
  const [loadingMessages, setLoadingMessages] = useState(false)
  const [sendingMessage, setSendingMessage] = useState(false)
  const [newModalOpen, setNewModalOpen] = useState(false)
  const [renameModalOpen, setRenameModalOpen] = useState(false)
  const [modelConfigOpen, setModelConfigOpen] = useState(false)
  const [modelConfigLoading, setModelConfigLoading] = useState(false)
  const [selectedKeyType, setSelectedKeyType] = useState<AgentApiKeyType>('OPENAI_API_KEY')
  const [apiKeyInput, setApiKeyInput] = useState('')
  const [newSessionName, setNewSessionName] = useState('New Bio Analysis')
  const [renameSessionName, setRenameSessionName] = useState('')
  const [messageApi, contextHolder] = message.useMessage()

  const loadSessions = useCallback(async () => {
    if (!projectId) return
    setLoadingSessions(true)
    try {
      const data = await listProjectSessions(projectId)
      if (data.length === 0) {
        const created = await createChatSession({
          project_id: projectId,
          initial_name: 'New Bio Analysis',
        })
        const fallbackSession: ChatSessionSummary = {
          session_id: created.session_id,
          name: created.name || 'New Bio Analysis',
        }
        setSessions([fallbackSession])
        setActiveSessionId(fallbackSession.session_id)
        return
      }
      setSessions(data)
      setActiveSessionId((prev) => (prev && data.some((s) => s.session_id === prev) ? prev : data[0].session_id))
    } catch (error) {
      const msg = error instanceof Error ? error.message : 'Failed to load session list'
      messageApi.error(msg)
      setSessions([])
      setActiveSessionId('')
    } finally {
      setLoadingSessions(false)
    }
  }, [projectId, messageApi])

  useEffect(() => {
    void loadSessions()
  }, [loadSessions])

  const loadMessages = useCallback(async (sessionId: string) => {
    if (!sessionId) {
      setMessages([])
      return
    }
    setLoadingMessages(true)
    try {
      const history = await getSessionHistory(sessionId, 1, 100)
      const mapped: Message[] = history.messages
        .map((item) => {
          const rawText = extractHistoryMessageText(item.content)
          const role: Message['role'] = item.type === 'user_input' ? 'user' : 'assistant'
          const text = role === 'assistant'
            ? normalizeAgentReply(rawText)
            : rawText.replace(/\\n/g, '\n')
          if (!text) return null
          return {
            id: item.id,
            role,
            content: text,
          }
        })
        .filter((x): x is Message => Boolean(x))
      setMessages([...mapped].reverse())
    } catch (error) {
      const msg = error instanceof Error ? error.message : 'Failed to load messages'
      messageApi.error(msg)
      setMessages([])
    } finally {
      setLoadingMessages(false)
    }
  }, [messageApi])

  useEffect(() => {
    void loadMessages(activeSessionId)
  }, [activeSessionId, loadMessages])

  const handleCreateSession = async (initialName: string) => {
    if (!projectId) return
    setCreatingSession(true)
    try {
      const created = await createChatSession({
        project_id: projectId,
        initial_name: initialName,
      })
      await loadSessions()
      if (created.session_id) {
        setActiveSessionId(created.session_id)
      }
      setMessages([])
    } catch (error) {
      const msg = error instanceof Error ? error.message : 'Failed to create session'
      messageApi.error(msg)
      throw error
    } finally {
      setCreatingSession(false)
    }
  }

  const handleRenameSession = async (sessionId: string, newName: string) => {
    if (!sessionId) return
    setRenamingSession(true)
    
    try {
      await renameChatSession(sessionId, { new_name: newName })
      await loadSessions()
    } catch (error) {
      const msg = error instanceof Error ? error.message : 'Failed to rename session'
      messageApi.error(msg)
      throw error
    } finally {
      setRenamingSession(false)
    }
  }

  const handleDeleteSession = async (sessionId: string) => {
    if (!sessionId) return
    setDeletingSession(true)
    try {
      await deleteChatSession(sessionId)
      await loadSessions()
      setMessages([])
    } catch (error) {
      const msg = error instanceof Error ? error.message : 'Failed to delete session'
      messageApi.error(msg)
      throw error
    } finally {
      setDeletingSession(false)
    }
  }

  const openRenameModal = () => {
    const active = sessions.find((s) => s.session_id === activeSessionId)
    if (!active) {
      messageApi.warning('No active session to rename.')
      return
    }
    setRenameSessionName(active.name || '')
    setRenameModalOpen(true)
  }

  const openDeleteConfirm = () => {
    if (!activeSessionId) {
      messageApi.warning('No active session to delete.')
      return
    }
    void Modal.confirm({
      title: 'Delete this session?',
      content: 'This action cannot be undone.',
      okText: 'Delete',
      cancelText: 'Cancel',
      okButtonProps: { danger: true },
      onOk: async () => {
        await handleDeleteSession(activeSessionId)
      },
    })
  }

  const handleSend = () => {
    if (!inputValue.trim()) return
    if (!activeSessionId) {
      messageApi.warning('No active session. Please create a session first.')
      return
    }
    const text = inputValue.trim()
    const userMessageId = `local-user-${Date.now()}`
    const thinkingMessageId = `local-thinking-${Date.now()}`
    setInputValue('')
    setMessages((prev) => ([
      ...prev,
      { id: userMessageId, role: 'user', content: text },
      { id: thinkingMessageId, role: 'assistant', content: 'Thinking...', pending: true },
    ]))
    void (async () => {
      setSendingMessage(true)
      try {
        const response = await agentChat({
          message: text,
          snapshot_id: snapshotId ?? '',
          project_id: projectId ?? '',
          dataset_id: datasetId ?? null,
          session_id: activeSessionId || null,
        })
        const displayReply = normalizeAgentReply(response.reply || '')
        setMessages((prev) => prev.map((msg) => (
          msg.id === thinkingMessageId
            ? { ...msg, content: displayReply || '(No response)', pending: false }
            : msg
        )))
        onAgentResponse?.()
      } catch (error) {
        const msg = error instanceof Error ? error.message : 'Failed to send message'
        setMessages((prev) => prev.map((item) => (
          item.id === thinkingMessageId
            ? { ...item, content: `Error: ${msg}`, pending: false }
            : item
        )))
        messageApi.error(msg)
      } finally {
        setSendingMessage(false)
      }
    })()
  }

  const handleSaveModelConfig = async () => {
    const keyValue = apiKeyInput.trim()
    if (!keyValue) {
      messageApi.warning('Please enter an API key')
      return
    }
    setModelConfigLoading(true)
    try {
      const verifyResult = await verifyApiKey({
        key_type: selectedKeyType,
        key_value: keyValue,
      })
      if (!verifyResult.valid) {
        throw new Error(verifyResult.message || 'Save failed')
      }

      await setApiKey({
        key_type: selectedKeyType,
        key_value: keyValue,
        is_current: true,
      })

      messageApi.success('Model config saved')
      setModelConfigOpen(false)
      setApiKeyInput('')
    } catch (error) {
      const msg = error instanceof Error ? `Save failed: ${error.message}` : 'Save failed'
      messageApi.error(msg)
    } finally {
      setModelConfigLoading(false)
    }
  }

  const handleKeyDown = (e: React.KeyboardEvent) => {
    if (e.key === 'Enter' && !e.shiftKey) {
      e.preventDefault()
      handleSend()
    }
  }

  return (
    <div className={`chat-container ${showHeader ? 'with-border' : ''}`}>
      {contextHolder}
      {showHeader && (
        <div className="chat-header">
          <Select
            value={activeSessionId || undefined}
            onChange={setActiveSessionId}
            placeholder="Select session"
            options={sessions.map((session) => ({
              value: session.session_id,
              label: session.name,
            }))}
            loading={loadingSessions}
            className="chat-header-session-select"
          />
          <div className="chat-header-actions">
            {showHeader && (
              <Dropdown
                trigger={['click']}
                menu={{
                  items: [
                    { key: 'model_config', label: 'Model Config' },
                    { key: 'new', label: 'New Session' },
                    { key: 'rename', label: 'Rename Session', disabled: !activeSessionId },
                    { key: 'delete', label: 'Delete Session', disabled: !activeSessionId, danger: true },
                  ],
                  onClick: ({ key }) => {
                    if (key === 'new') {
                      setNewSessionName('New Bio Analysis')
                      setNewModalOpen(true)
                    } else if (key === 'rename') {
                      openRenameModal()
                    } else if (key === 'model_config') {
                      setModelConfigOpen(true)
                    } else if (key === 'delete') {
                      openDeleteConfirm()
                    }
                  },
                }}
              >
                <Button
                  type="text"
                  shape="circle"
                  icon={<EllipsisOutlined />}
                  className="chat-ellipsis-btn"
                  title="Session Actions"
                />
              </Dropdown>
            )}
            {onClose && (
              <Button
                type="text"
                icon={<CloseOutlined />}
                onClick={onClose}
                className="chat-close-btn"
                title="Hide Chat"
              />
            )}
          </div>
        </div>
      )}
      {!showHeader && (
        <SessionToolbar
          sessions={sessions}
          activeSessionId={activeSessionId}
          loadingSessions={loadingSessions}
          creatingSession={creatingSession}
          renamingSession={renamingSession}
          deletingSession={deletingSession}
          hideActions={false}
          onOpenModelConfig={() => setModelConfigOpen(true)}
          onSelectSession={setActiveSessionId}
          onCreateSession={handleCreateSession}
          onRenameSession={handleRenameSession}
          onDeleteSession={handleDeleteSession}
        />
      )}

      <div className="messages-area">
        {messages.length === 0 ? (
          <div className="empty-state">
            Start a conversation
          </div>
        ) : (
          messages.map(msg => (
            <div key={msg.id} className={`message ${msg.role}${msg.pending ? ' thinking' : ''}`}>
              {msg.role === 'assistant'
                ? <ChatMessageContent content={msg.content} />
                : msg.content}
            </div>
          ))
        )}
      </div>

      <div className="input-area">
        <div className="input-wrapper">
          <TextArea
            value={inputValue}
            onChange={e => setInputValue(e.target.value)}
            onKeyDown={handleKeyDown}
            placeholder="Message..."
            autoSize={{ minRows: 1, maxRows: 6 }}
            className="message-input"
          />
          <Button
            type="primary"
            icon={<SendOutlined />}
            onClick={handleSend}
            loading={sendingMessage}
            disabled={loadingMessages || !activeSessionId || !projectId}
            className="send-btn"
          />
        </div>
      </div>
       <Modal
        title="Model Config"
        open={modelConfigOpen}
        onCancel={() => setModelConfigOpen(false)}
        onOk={handleSaveModelConfig}
        okText="Save"
        confirmLoading={modelConfigLoading}
      >
        <div style={{ display: 'flex', flexDirection: 'column', gap: 12 }}>
          <div>Model</div>
          <Select
            value={selectedKeyType}
            onChange={(value) => setSelectedKeyType(value as AgentApiKeyType)}
            showSearch
            optionFilterProp="label"
            options={MODEL_KEY_OPTIONS.map((option) => ({
              value: option,
              label: getModelKeyDisplayLabel(option),
            }))}
            placeholder="Select model key type"
          />
          <div>API KEY</div>
          <AntInput.Password
            value={apiKeyInput}
            onChange={(e) => setApiKeyInput(e.target.value)}
            placeholder="Enter API key"
          />
        </div>
      </Modal>
      <Modal
        title="Create New Session"
        open={newModalOpen}
        onCancel={() => setNewModalOpen(false)}
        onOk={async () => {
          const next = newSessionName.trim()
          if (!next) return
          await handleCreateSession(next)
          setNewModalOpen(false)
        }}
        okText="Create"
        confirmLoading={creatingSession}
      >
        <AntInput
          value={newSessionName}
          onChange={(e) => setNewSessionName(e.target.value)}
          placeholder="Session name"
          maxLength={120}
        />
      </Modal>

      <Modal
        title="Rename Session"
        open={renameModalOpen}
        onCancel={() => setRenameModalOpen(false)}
        onOk={async () => {
          const next = renameSessionName.trim()
          if (!next || !activeSessionId) return
          await handleRenameSession(activeSessionId, next)
          setRenameModalOpen(false)
        }}
        okText="Save"
        confirmLoading={renamingSession}
      >
        <AntInput
          value={renameSessionName}
          onChange={(e) => setRenameSessionName(e.target.value)}
          placeholder="Session name"
          maxLength={120}
        />
      </Modal>

    </div>
  )
}
