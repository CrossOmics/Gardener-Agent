import { useEffect, useMemo, useState } from 'react'
import { Alert, Button, Empty, Spin, Table } from 'antd'
import { getSnapshotAncestors } from '../../../services/endpoints'
import type { SnapshotAncestorNode } from '../../../services/types'

interface LineagePanelProps {
  snapshotId: string
  onNavigateToSnapshot?: (snapshotId: string) => void
}

interface LineageRow {
  key: string
  snapshotId: string
  snapshotName: string
  stageName: string
  createTime: string
}

function formatCreateTime(raw: string): string {
  const text = (raw || '').trim()
  if (!text) return 'N/A'
  const match = text.match(/^(\d{4})-(\d{2})-(\d{2})T(\d{2}):(\d{2}):(\d{2})/)
  if (!match) {
    return text.split('.')[0].replace('T', ' ')
  }

  const [, y, m, d, hh, mm, ss] = match
  const utcMs = Date.UTC(
    Number(y),
    Number(m) - 1,
    Number(d),
    Number(hh),
    Number(mm),
    Number(ss)
  )
  const offsetMinutes = new Date().getTimezoneOffset()
  const shifted = new Date(utcMs - offsetMinutes * 60 * 1000)

  const pad = (value: number) => String(value).padStart(2, '0')
  const year = shifted.getUTCFullYear()
  const month = pad(shifted.getUTCMonth() + 1)
  const day = pad(shifted.getUTCDate())
  const hour = pad(shifted.getUTCHours())
  const minute = pad(shifted.getUTCMinutes())
  const second = pad(shifted.getUTCSeconds())
  return `${year}-${month}-${day} ${hour}:${minute}:${second}`
}

export default function LineagePanel({
  snapshotId,
  onNavigateToSnapshot,
}: LineagePanelProps) {
  const [loading, setLoading] = useState(false)
  const [error, setError] = useState<string | null>(null)
  const [ancestors, setAncestors] = useState<SnapshotAncestorNode[]>([])

  useEffect(() => {
    let cancelled = false
    const load = async () => {
      if (!snapshotId) {
        setAncestors([])
        setError(null)
        return
      }
      setLoading(true)
      setError(null)
      try {
        const data = await getSnapshotAncestors(snapshotId)
        if (!cancelled) {
          setAncestors(Array.isArray(data) ? data : [])
        }
      } catch (e) {
        if (!cancelled) {
          setAncestors([])
          setError(e instanceof Error ? e.message : 'Failed to load lineage')
        }
      } finally {
        if (!cancelled) {
          setLoading(false)
        }
      }
    }
    void load()
    return () => {
      cancelled = true
    }
  }, [snapshotId])

  const rows = useMemo<LineageRow[]>(
    () =>
      [...ancestors].reverse().map((node, index) => ({
        key: `${node.snapshot_id}-${index}`,
        snapshotId: node.snapshot_id,
        snapshotName: node.snapshot_name && node.snapshot_name.trim() ? node.snapshot_name : 'Null',
        stageName: node.stage_name || 'Unknown',
        createTime: formatCreateTime(node.create_time || ''),
      })),
    [ancestors]
  )

  return (
    <div className="lineage-panel">
      {error && <Alert type="error" showIcon message={error} />}
      {loading ? (
        <div className="lineage-loading">
          <Spin size="small" />
        </div>
      ) : rows.length === 0 ? (
        <Empty image={Empty.PRESENTED_IMAGE_SIMPLE} description="No lineage data" />
      ) : (
        <Table<LineageRow>
          size="small"
          pagination={false}
          rowKey="key"
          dataSource={rows}
          columns={[
            {
              title: 'Snapshot Name',
              dataIndex: 'snapshotName',
              key: 'snapshotName',
              width: 220,
            },
            {
              title: 'Stage',
              dataIndex: 'stageName',
              key: 'stageName',
              width: 140,
            },
            {
              title: 'Created At',
              dataIndex: 'createTime',
              key: 'createTime',
              width: 220,
            },
            {
              title: 'Jump to the folder',
              key: 'jump',
              width: 92,
              render: (_, record) => (
                <Button
                  size="small"
                  onClick={() => onNavigateToSnapshot?.(record.snapshotId)}
                >
                  Go
                </Button>
              ),
            },
          ]}
        />
      )}
    </div>
  )
}
