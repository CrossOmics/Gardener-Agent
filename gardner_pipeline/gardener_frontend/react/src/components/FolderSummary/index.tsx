import { useCallback, useEffect, useMemo, useState } from 'react'
import { Typography, Divider, Form, Input, Switch, Spin } from 'antd'
import { FolderOpenOutlined } from '@ant-design/icons'
import { getCategoryFromFolderName, isPreprocessingSubfolder } from './folderCategory'
import { getProjectDisplayMap, getSnapshotDetails } from '../../services/endpoints'
import FolderSettingsPanel from './FolderSettingsPanel'
import './styles.css'

const { Title, Text } = Typography

interface FolderSummaryProps {
  folderPath: string
  folderName: string
  projectId?: string
  datasetId?: string
  snapshotId?: string
  isSnapshotFolder?: boolean
  snapshotCategory?: 'preprocessing' | 'clustering' | 'dge' | 'annotation'
  onRunSuccess?: () => void
}

export default function FolderSummary({
  folderPath,
  folderName,
  projectId,
  datasetId,
  snapshotId,
  isSnapshotFolder: isSnapshotFolderProp,
  snapshotCategory: snapshotCategoryProp,
  onRunSuccess,
}: FolderSummaryProps) {
  const category = getCategoryFromFolderName(folderName)
  const [snapshotLoading, setSnapshotLoading] = useState(false)
  const [snapshotError, setSnapshotError] = useState<string | null>(null)
  const [snapshotInput, setSnapshotInput] = useState<Record<string, unknown> | null>(null)
  const [parentSnapshotId, setParentSnapshotId] = useState<string | null>(null)
  const [parentSnapshotName, setParentSnapshotName] = useState<string | null>(null)

  const findSnapshotDisplayName = useCallback((id: string, snapshotsByStage: Record<string, Record<string, string>>) => {
    for (const stageMap of Object.values(snapshotsByStage)) {
      if (stageMap[id]) return stageMap[id]
    }
    return null
  }, [])

  const isSnapshotFolder = Boolean(isSnapshotFolderProp)
  // Don't show editable settings panel for snapshot folders or preprocessing subfolders.
  const shouldShowSettings = category && projectId && datasetId && !isPreprocessingSubfolder(folderName) && !isSnapshotFolder

  const snapshotCategory = snapshotCategoryProp ?? null

  useEffect(() => {
    let cancelled = false
    const loadSnapshot = async () => {
      if (!isSnapshotFolder) {
        setSnapshotInput(null)
        setParentSnapshotId(null)
        setParentSnapshotName(null)
        setSnapshotError(null)
        return
      }
      setSnapshotLoading(true)
      setSnapshotError(null)
      try {
        const detail = await getSnapshotDetails(snapshotId || folderName)
        const input = detail.input && typeof detail.input === 'object' ? (detail.input as Record<string, unknown>) : {}
        if (!cancelled) {
          setSnapshotInput(input)
          const parentId = detail.parent_snapshot_id ?? null
          setParentSnapshotId(parentId)
          if (parentId && projectId && datasetId) {
            try {
              const displayMap = await getProjectDisplayMap(projectId)
              const datasetDisplay = displayMap.datasets?.[datasetId]
              const displayName = datasetDisplay ? findSnapshotDisplayName(parentId, datasetDisplay.snapshots) : null
              setParentSnapshotName(displayName ?? null)
            } catch {
              setParentSnapshotName(null)
            }
          } else {
            setParentSnapshotName(null)
          }
        }
      } catch (error) {
        if (!cancelled) {
          setSnapshotInput(null)
          setParentSnapshotId(null)
          setParentSnapshotName(null)
          setSnapshotError(error instanceof Error ? error.message : 'Failed to load snapshot details')
        }
      } finally {
        if (!cancelled) {
          setSnapshotLoading(false)
        }
      }
    }
    void loadSnapshot()
    return () => {
      cancelled = true
    }
  }, [isSnapshotFolder, snapshotId, folderName, projectId, datasetId, findSnapshotDisplayName])

  const selectedMethods = useMemo(() => {
    if (!snapshotInput) return [] as Array<{ id: string; name: string; type: 'celltypist' | 'gseapy' }>
    if (snapshotCategory !== 'annotation') return [] as Array<{ id: string; name: string; type: 'celltypist' | 'gseapy' }>
    const modelNames = Array.isArray(snapshotInput.model_names) ? snapshotInput.model_names.map((v) => String(v)) : []
    const categories = Array.isArray(snapshotInput.categories) ? snapshotInput.categories.map((v) => String(v)) : []
    return [
      ...modelNames.map((name) => ({ id: `ct:${name}`, name, type: 'celltypist' as const })),
      ...categories.map((name) => ({ id: `gs:${name}`, name, type: 'gseapy' as const })),
    ]
  }, [snapshotInput, snapshotCategory])

  const snapshotReadonlyValues = useMemo(() => {
    if (!snapshotInput || !snapshotCategory) return null
    if (snapshotCategory === 'preprocessing') {
      return {
        // keep tab-based readonly rendering, but align values with snapshots/query input
        skip_qc_calculation: snapshotInput.skip_qc_calculation ?? false,
        skip_qc_filter: snapshotInput.skip_qc_filter ?? false,
        min_genes: snapshotInput.min_genes,
        min_cells: snapshotInput.min_cells,
        pct_mt_max: snapshotInput.pct_mt_max ?? null,
        max_counts: snapshotInput.max_counts ?? snapshotInput.cell_max_counts ?? null,
        pct_hb_max: snapshotInput.pct_hb_max ?? null,
        skip_hvg: snapshotInput.skip_hvg ?? false,
        n_top_genes_hvg: snapshotInput.n_top_genes,
        n_top_genes: snapshotInput.n_top_genes,
        flavor: snapshotInput.flavor,
        target_sum: snapshotInput.target_sum,
        skip_pca: snapshotInput.skip_pca ?? false,
        n_comps: snapshotInput.n_comps,
        svd_solver: snapshotInput.svd_solver,
        skip_neighbors: snapshotInput.skip_neighbors ?? false,
        n_neighbors: snapshotInput.n_neighbors,
        n_pcs: snapshotInput.n_pcs,
      }
    }
    if (snapshotCategory === 'clustering') {
      return {
        snapshot_id: parentSnapshotName ?? 'Null',
        method: snapshotInput.method,
        resolution: snapshotInput.resolution,
        run_hierarchical: snapshotInput.run_hierarchical,
      }
    }
    if (snapshotCategory === 'dge') {
      return {
        snapshot_id: parentSnapshotName ?? 'Null',
        groupby: snapshotInput.groupby,
        method: snapshotInput.method,
        n_top_genes: snapshotInput.n_top_genes,
        use_raw: snapshotInput.use_raw,
      }
    }
    return {
      snapshot_id: parentSnapshotName ?? 'Null',
      annotation_majority_voting: snapshotInput.majority_voting ?? true,
      annotation_top_n_genes: snapshotInput.top_n_genes ?? 100,
    }
  }, [snapshotInput, snapshotCategory, parentSnapshotId, parentSnapshotName])

  return (
    <div className="folder-summary-container">
      <div className="folder-summary-header">
        <FolderOpenOutlined className="folder-summary-icon" />
        <Title level={3} className="folder-summary-title">{folderName}</Title>
      </div>

      {/* Settings panel for matching category folders */}
      {shouldShowSettings && (
        <>
          <FolderSettingsPanel
            projectId={projectId}
            datasetId={datasetId}
            folderPath={folderPath}
            folderName={folderName}
            category={category}
            onRunSuccess={onRunSuccess}
          />
        </>
      )}

      {isSnapshotFolder && (
        <>
          <Divider />
          {snapshotLoading ? (
            <div className="folder-settings-form snapshot-readonly">
              <div className="snapshot-loading-wrap">
                <Spin size="small" />
              </div>
            </div>
          ) : snapshotError ? (
            <div className="folder-settings-form snapshot-readonly">
              <Text type="secondary">{snapshotError}</Text>
            </div>
          ) : snapshotCategory && snapshotReadonlyValues ? (
            <FolderSettingsPanel
              projectId={projectId || ''}
              datasetId={datasetId || ''}
              folderPath={folderPath}
              folderName={folderName}
              category={snapshotCategory}
              readOnly
              readOnlyValues={snapshotReadonlyValues}
              readOnlySelectedMethods={selectedMethods}
            />
          ) : (
            <div className="folder-settings-form snapshot-readonly">
              <Form layout="vertical" size="small">
                {Object.entries(snapshotInput || {}).map(([key, value]) => (
                  <Form.Item key={key} label={key}>
                    {typeof value === 'boolean' ? (
                      <Switch checked={value} disabled />
                    ) : typeof value === 'number' || typeof value === 'string' ? (
                      <Input value={String(value)} disabled />
                    ) : (
                      <Input.TextArea
                        value={JSON.stringify(value, null, 2)}
                        autoSize={{ minRows: 2, maxRows: 6 }}
                        disabled
                      />
                    )}
                  </Form.Item>
                ))}
                {snapshotCategory === 'annotation' && (
                  <div className="selected-section">
                    <div className="selected-header">
                      <h3 className="selected-title">Selected Methods</h3>
                      <span className="color-legend">
                        <span className="legend-item">
                          <span className="legend-dot celltypist-dot"></span>CellTypist
                        </span>
                        <span className="legend-item">
                          <span className="legend-dot gseapy-dot"></span>GSEApy
                        </span>
                      </span>
                    </div>
                    <div className="selected-models">
                      {selectedMethods.length === 0 ? (
                        <span className="no-selection">No methods selected</span>
                      ) : (
                        selectedMethods.map((model) => (
                          <div
                            key={model.id}
                            className="selected-tag"
                            style={{ backgroundColor: model.type === 'celltypist' ? '#fa8c16' : '#1890ff' }}
                          >
                            <span className="tag-name">{model.name}</span>
                          </div>
                        ))
                      )}
                    </div>
                  </div>
                )}
              </Form>
            </div>
          )}
        </>
      )}

    </div>
  )
}
