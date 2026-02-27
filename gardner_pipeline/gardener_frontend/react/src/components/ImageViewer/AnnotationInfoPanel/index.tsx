import { useEffect, useMemo, useState } from 'react'
import { Button, Input, Modal, Space, Table, Typography, message } from 'antd'
import { getSnapshotDetails, updateAnnotationLabels } from '../../../services/endpoints'
import './styles.css'

const { Text } = Typography

interface AnnotationInfoPanelProps {
  imagePath: string
  snapshotId: string
  parentSnapshotId?: string | null
  snapshotOutput: Record<string, unknown> | null
  onNavigateToSnapshot?: (snapshotId: string) => void
  onRenameSuccess?: () => void
}

interface AnnotationRow {
  id: string
  predictedCellType: string
  averageConfidence: number | null
  color: string
}

function scopedNameKey(annotationKey: string | null, clusterId: string): string {
  return `${annotationKey ?? '__unknown__'}::${clusterId}`
}

function toImageKey(path: string): string {
  const normalized = path.replace(/\\/g, '/')
  const fileName = normalized.split('/').pop() || ''
  return fileName.replace(/\.[^/.]+$/, '')
}

function toRow(id: string, value: unknown): AnnotationRow | null {
  if (!value || typeof value !== 'object') {
    return null
  }
  const record = value as Record<string, unknown>
  const predictedCellType = typeof record.predicted_cell_type === 'string' ? record.predicted_cell_type : ''
  const averageConfidence = typeof record.average_confidence === 'number' ? record.average_confidence : null
  const color = typeof record.color === 'string' ? record.color : '#d9d9d9'
  return { id, predictedCellType, averageConfidence, color }
}

export default function AnnotationInfoPanel({
  imagePath,
  snapshotId,
  parentSnapshotId,
  snapshotOutput,
  onNavigateToSnapshot,
  onRenameSuccess,
}: AnnotationInfoPanelProps) {
  const [messageApi, contextHolder] = message.useMessage()
  const [editableNames, setEditableNames] = useState<Record<string, string>>({})
  const [renameLoading, setRenameLoading] = useState(false)
  const [detailOpen, setDetailOpen] = useState(false)
  const [detailLoading, setDetailLoading] = useState(false)
  const [detailJson, setDetailJson] = useState('{}')
  const [pagination, setPagination] = useState({ current: 1, pageSize: 5 })

  const imageKey = useMemo(() => toImageKey(imagePath), [imagePath])

  const matchedOutputKey = useMemo(() => {
    if (!snapshotOutput) return null
    const keys = Object.keys(snapshotOutput)
    const exact = keys.find((key) => key === imageKey)
    if (exact) return exact
    const lower = imageKey.toLowerCase()
    return keys.find((key) => key.toLowerCase() === lower) ?? null
  }, [imageKey, snapshotOutput])

  const rows = useMemo(() => {
    if (!matchedOutputKey || !snapshotOutput) return [] as AnnotationRow[]
    const target = snapshotOutput[matchedOutputKey]
    if (!target || typeof target !== 'object') return [] as AnnotationRow[]
    const mapped = Object.entries(target as Record<string, unknown>)
      .map(([id, value]) => toRow(id, value))
      .filter((row): row is AnnotationRow => row !== null)

    return mapped.sort((a, b) => {
      const aNum = Number(a.id)
      const bNum = Number(b.id)
      if (Number.isNaN(aNum) || Number.isNaN(bNum)) return a.id.localeCompare(b.id)
      return aNum - bNum
    })
  }, [matchedOutputKey, snapshotOutput])

  useEffect(() => {
    setEditableNames((prev) => {
      if (rows.length === 0) {
        return prev
      }
      const next: Record<string, string> = {}
      for (const row of rows) {
        const key = scopedNameKey(matchedOutputKey, row.id)
        next[key] = prev[key] ?? row.predictedCellType
      }
      const prevKeys = Object.keys(prev).filter((key) => key.startsWith(`${matchedOutputKey ?? '__unknown__'}::`))
      const nextKeys = Object.keys(next)
      if (prevKeys.length === nextKeys.length && nextKeys.every((key) => prev[key] === next[key])) {
        return prev
      }
      return { ...prev, ...next }
    })
  }, [rows, matchedOutputKey])

  useEffect(() => {
    setPagination((prev) => ({ ...prev, current: 1 }))
  }, [matchedOutputKey, imagePath])

  const handleOpenDetail = async (clusterId: string) => {
    if (!parentSnapshotId) {
      messageApi.error('No parent DGE snapshot found')
      return
    }

    setDetailOpen(true)
    setDetailLoading(true)
    setDetailJson('{}')

    try {
      const dgeDetail = await getSnapshotDetails(parentSnapshotId)
      const output = dgeDetail.output && typeof dgeDetail.output === 'object'
        ? (dgeDetail.output as Record<string, unknown>)
        : {}
      const topMarkers = output.top_markers && typeof output.top_markers === 'object'
        ? (output.top_markers as Record<string, unknown>)
        : {}
      const focused = topMarkers[clusterId]
      const payload = focused !== undefined ? { [clusterId]: focused } : topMarkers
      setDetailJson(JSON.stringify(payload, null, 2))
    } catch (error) {
      const msg = error instanceof Error ? error.message : 'Failed to load DGE markers'
      setDetailJson(JSON.stringify({ error: msg }, null, 2))
    } finally {
      setDetailLoading(false)
    }
  }

  const handleRename = async () => {
    if (!matchedOutputKey) {
      messageApi.error('No matched annotation key for current image')
      return
    }

    const changed: Record<string, string> = {}
    for (const row of rows) {
      const value = (editableNames[scopedNameKey(matchedOutputKey, row.id)] ?? '').trim()
      if (value && value !== row.predictedCellType) {
        changed[row.id] = value
      }
    }

    if (Object.keys(changed).length === 0) {
      messageApi.info('No changed labels to update')
      return
    }

    setRenameLoading(true)
    try {
      const resp = await updateAnnotationLabels({
        snapshot_id: snapshotId,
        file_name: matchedOutputKey,
        updated_annotation: changed,
      })
      messageApi.success(resp.msg || 'Annotation labels updated')
      onRenameSuccess?.()
    } catch (error) {
      messageApi.error(error instanceof Error ? error.message : 'Update annotation labels failed')
    } finally {
      setRenameLoading(false)
    }
  }

  const handleResetNames = () => {
    const reset: Record<string, string> = { ...editableNames }
    for (const row of rows) {
      reset[scopedNameKey(matchedOutputKey, row.id)] = row.predictedCellType
    }
    setEditableNames(reset)
    messageApi.success('Names reset to original values')
  }

  return (
    <div className="annotation-info-panel">
      {contextHolder}
      {!matchedOutputKey ? (
        <Text type="secondary">No matching annotation output for image key: {imageKey}</Text>
      ) : (
        <>
          <Table<AnnotationRow>
            rowKey="id"
            size="small"
            pagination={{
              current: pagination.current,
              pageSize: pagination.pageSize,
              showSizeChanger: true,
              pageSizeOptions: [5, 10, 20],
              size: 'small',
              onChange: (current, pageSize) => {
                setPagination({
                  current,
                  pageSize: pageSize || 5,
                })
              },
            }}
            dataSource={rows}
            className="annotation-info-table"
            columns={[
              {
                title: 'ID',
                dataIndex: 'id',
                width: 80,
              },
              {
                title: 'Predicted Cell Type',
                dataIndex: 'predictedCellType',
                render: (_, record) => (
                  <Input
                    value={editableNames[scopedNameKey(matchedOutputKey, record.id)] ?? record.predictedCellType}
                    onChange={(e) => {
                      const next = e.target.value
                      setEditableNames((prev) => ({
                        ...prev,
                        [scopedNameKey(matchedOutputKey, record.id)]: next,
                      }))
                    }}
                  />
                ),
              },
              {
                title: 'Average Confidence',
                dataIndex: 'averageConfidence',
                width: 180,
                render: (value: number | null) => (value === null ? 'N/A' : value.toFixed(4)),
              },
              {
                title: 'Color',
                dataIndex: 'color',
                width: 120,
                render: (value: string) => <span className="annotation-color-dot" style={{ backgroundColor: value }} />,
              },
              {
                title: 'Detail',
                key: 'detail',
                width: 120,
                render: (_, record) => (
                  <Button type="link" onClick={() => void handleOpenDetail(record.id)}>
                    Detail
                  </Button>
                ),
              },
            ]}
          />
          <div className="annotation-info-actions">
            <Button onClick={handleResetNames} disabled={rows.length === 0}>
              Reset Names
            </Button>
            <Button type="primary" onClick={() => void handleRename()} loading={renameLoading}>
              Rename
            </Button>
          </div>
        </>
      )}

      <Modal
        title="DGE Top Markers"
        open={detailOpen}
        onCancel={() => setDetailOpen(false)}
        footer={[
          <Button key="close" onClick={() => setDetailOpen(false)}>
            Close
          </Button>,
          <Button
            key="go"
            type="primary"
            disabled={!parentSnapshotId || !onNavigateToSnapshot}
            onClick={() => {
              if (!parentSnapshotId || !onNavigateToSnapshot) return
              onNavigateToSnapshot(parentSnapshotId)
              setDetailOpen(false)
            }}
          >
            Go to DGE Snapshot
          </Button>,
        ]}
      >
        {detailLoading ? (
          <Space direction="vertical">
            <Text type="secondary">Loading...</Text>
          </Space>
        ) : (
          <pre className="annotation-detail-pre">{detailJson}</pre>
        )}
      </Modal>
    </div>
  )
}
