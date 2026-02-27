import { useMemo, useState } from 'react'
import { Button, Checkbox, Input, message, Switch, Tooltip } from 'antd'
import { InfoCircleOutlined } from '@ant-design/icons'
import { runClusteringMerge, runClusteringSubcluster } from '../../../services/endpoints'
import './styles.css'

interface ClusterActionsTableProps {
  projectId?: string
  datasetId?: string
  snapshotId?: string | null
  method?: 'leiden' | 'louvain' | 'cplearn'
  resolution?: number
  clusterIds: Record<string, string>
  onRunSuccess?: () => void
}

export default function ClusterActionsTable({
  projectId,
  datasetId,
  snapshotId,
  method = 'leiden',
  resolution,
  clusterIds,
  onRunSuccess,
}: ClusterActionsTableProps) {
  const [selected, setSelected] = useState<Set<string>>(new Set())
  const [mergeLabel, setMergeLabel] = useState('')
  const [overwrite, setOverwrite] = useState(false)
  const [loadingAction, setLoadingAction] = useState<'merge' | 'cluster' | null>(null)
  const [messageApi, contextHolder] = message.useMessage()

  const orderedEntries = useMemo(() => {
    return Object.entries(clusterIds).sort(([a], [b]) => {
      const an = Number(a)
      const bn = Number(b)
      if (Number.isFinite(an) && Number.isFinite(bn)) return an - bn
      return a.localeCompare(b)
    })
  }, [clusterIds])

  const toggle = (id: string) => {
    setSelected((prev) => {
      const next = new Set(prev)
      if (next.has(id)) next.delete(id)
      else next.add(id)
      return next
    })
  }

  const ensureContext = () => {
    if (!projectId || !datasetId || !snapshotId) {
      messageApi.error('Missing project/dataset/snapshot context')
      return false
    }
    return true
  }

  const handleMerge = async () => {
    if (!ensureContext()) return
    const ids = Array.from(selected)
    if (ids.length < 2) {
      messageApi.warning('Select at least two clusters for Merge')
      return
    }
    setLoadingAction('merge')
    try {
      const requestBody = {
        project_id: projectId!,
        dataset_id: datasetId!,
        snapshot_id: snapshotId!,
        method,
        clusters_to_merge: ids,
        overwrite,
        ...(mergeLabel.trim() ? { new_label: mergeLabel.trim() } : {}),
      }
      const response = await runClusteringMerge(requestBody)
      messageApi.success(response.msg || 'Merge completed')
      setMergeLabel('')
      onRunSuccess?.()
    } catch (error) {
      messageApi.error(error instanceof Error ? error.message : 'Merge failed')
    } finally {
      setLoadingAction(null)
    }
  }

  const handleContinueClustering = async () => {
    if (!ensureContext()) return
    const ids = Array.from(selected)
    if (ids.length === 0) {
      messageApi.warning('Select at least one cluster for Sub Clustering')
      return
    }
    setLoadingAction('cluster')
    try {
      const response = await runClusteringSubcluster({
        project_id: projectId!,
        dataset_id: datasetId!,
        snapshot_id: snapshotId!,
        source_cluster_col: method,
        target_clusters: ids,
        method,
        resolution,
        overwrite,
      })
      messageApi.success(response.msg || 'Sub-clustering completed')
      onRunSuccess?.()
    } catch (error) {
      messageApi.error(error instanceof Error ? error.message : 'Sub-clustering failed')
    } finally {
      setLoadingAction(null)
    }
  }

  if (orderedEntries.length === 0) {
    return (
      <div className="cluster-actions-empty">No cluster IDs available.</div>
    )
  }

  return (
    <div className="cluster-actions-wrap">
      {contextHolder}
      <div className="cluster-cards" aria-label="cluster ids cards">
        {orderedEntries.map(([id, color]) => (
          <div key={`card-${id}`} className="cluster-card">
            <div className="cluster-card-row">
              <span className="cluster-card-label">ID</span>
              <span className="cluster-card-value">{id}</span>
            </div>
            <div className="cluster-card-row">
              <span className="cluster-card-label">Color</span>
              <span className="cluster-color-dot" style={{ backgroundColor: color }} />
            </div>
            <div className="cluster-card-row">
              <span className="cluster-card-label">Selected</span>
              <Checkbox checked={selected.has(id)} onChange={() => toggle(id)} />
            </div>
          </div>
        ))}
      </div>

      <div className="cluster-actions-buttons">
        <div className="cluster-actions-params">
          <Input
            className="merge-label-input"
            placeholder="New cluster label (optional)"
            value={mergeLabel}
            onChange={(e) => setMergeLabel(e.target.value)}
            onPressEnter={() => {
              if (loadingAction !== 'merge') {
                void handleMerge()
              }
            }}
          />
          <div className="sub-overwrite-group">
            <span className="sub-overwrite-label">Overwrite</span>
            <Switch checked={overwrite} onChange={setOverwrite} />
            <Tooltip title={<span>On: overwrite current snapshot.<br />Off: create a new sub-clustering snapshot.</span>}>
              <InfoCircleOutlined className="sub-overwrite-tip" />
            </Tooltip>
          </div>
        </div>
        <div className="cluster-actions-row">
          <Button
            className="cluster-btn cluster-btn-merge"
            onClick={handleMerge}
            loading={loadingAction === 'merge'}
          >
            Merge
          </Button>
          <Button
            className="cluster-btn cluster-btn-sub"
            onClick={handleContinueClustering}
            loading={loadingAction === 'cluster'}
          >
            Sub Clustering
          </Button>
        </div>
      </div>
    </div>
  )
}
