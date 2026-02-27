import { useState } from 'react'
import { Button, Input, Modal, Popconfirm, Select } from 'antd'
import { PlusOutlined, EditOutlined, DeleteOutlined, SettingOutlined } from '@ant-design/icons'
import type { ChatSessionSummary } from '../../../services/types'
import './styles.css'

interface SessionToolbarProps {
  sessions: ChatSessionSummary[]
  activeSessionId: string
  loadingSessions: boolean
  creatingSession: boolean
  renamingSession: boolean
  deletingSession: boolean
  hideActions?: boolean
  onOpenModelConfig?: () => void
  onSelectSession: (sessionId: string) => void
  onCreateSession: (initialName: string) => Promise<void>
  onRenameSession: (sessionId: string, newName: string) => Promise<void>
  onDeleteSession: (sessionId: string) => Promise<void>
}

export default function SessionToolbar({
  sessions,
  activeSessionId,
  loadingSessions,
  creatingSession,
  renamingSession,
  deletingSession,
  hideActions = false,
  onOpenModelConfig,
  onSelectSession,
  onCreateSession,
  onRenameSession,
  onDeleteSession,
}: SessionToolbarProps) {
  const [newModalOpen, setNewModalOpen] = useState(false)
  const [renameModalOpen, setRenameModalOpen] = useState(false)
  const [newSessionName, setNewSessionName] = useState('New Bio Analysis')
  const [renameSessionName, setRenameSessionName] = useState('')

  const activeSession = sessions.find((session) => session.session_id === activeSessionId)

  const openRenameModal = () => {
    if (!activeSession) return
    setRenameSessionName(activeSession.name || '')
    setRenameModalOpen(true)
  }

  return (
    <div className="session-toolbar">
      <Select
        value={activeSessionId || undefined}
        onChange={onSelectSession}
        placeholder="Select session"
        options={sessions.map((session) => ({
          value: session.session_id,
          label: session.name,
        }))}
        loading={loadingSessions}
        className="session-select"
      />
      {!hideActions && (
        <>
          <Button
            icon={<PlusOutlined />}
            onClick={() => {
              setNewSessionName('New Bio Analysis')
              setNewModalOpen(true)
            }}
            loading={creatingSession}
            className="session-action-btn"
          >
            New
          </Button>
          <Button
            icon={<EditOutlined />}
            onClick={openRenameModal}
            disabled={!activeSessionId}
            loading={renamingSession}
            className="session-action-btn"
          >
            Rename
          </Button>
          <Popconfirm
            title="Delete this session?"
            description="This action cannot be undone."
            onConfirm={() => onDeleteSession(activeSessionId)}
            okText="Delete"
            cancelText="Cancel"
            disabled={!activeSessionId}
          >
            <Button
              icon={<DeleteOutlined />}
              disabled={!activeSessionId}
              loading={deletingSession}
              className="session-action-btn"
            >
              Delete
            </Button>
          </Popconfirm>
          <Button
            icon={<SettingOutlined />}
            onClick={onOpenModelConfig}
            className="session-action-btn"
          >
            Model Config
          </Button>
        </>
      )}

      {!hideActions && (
        <>
          <Modal
            title="Create New Session"
            open={newModalOpen}
            onCancel={() => setNewModalOpen(false)}
            onOk={async () => {
              const next = newSessionName.trim()
              if (!next) return
              await onCreateSession(next)
              setNewModalOpen(false)
            }}
            okText="Create"
            confirmLoading={creatingSession}
          >
            <Input
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
              await onRenameSession(activeSessionId, next)
              setRenameModalOpen(false)
            }}
            okText="Save"
            confirmLoading={renamingSession}
          >
            <Input
              value={renameSessionName}
              onChange={(e) => setRenameSessionName(e.target.value)}
              placeholder="Session name"
              maxLength={120}
            />
          </Modal>
        </>
      )}
    </div>
  )
}
