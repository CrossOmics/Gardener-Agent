import React, { useEffect, useRef, useState } from 'react'
import { Empty, Typography, Spin, Tabs } from 'antd'
import { FileImageOutlined } from '@ant-design/icons'
import { getSnapshotDetails } from '../../services/endpoints'
import { buildClusteringSummary, buildPreprocessingSummary } from './snapshotSummary'
import ClusterActionsTable from './ClusterActionsTable'
import AnnotationInfoPanel from './AnnotationInfoPanel'
import LineagePanel from './LineagePanel'
import './styles.css'

const { Text } = Typography

interface ImageViewerProps {
  imagePath: string | null
  snapshotId?: string | null
  projectId?: string
  datasetId?: string
  onZoomChange?: (zoom: number) => void
  onNavigateToSnapshot?: (snapshotId: string) => void
  onStructureChanged?: () => void
}

export default function ImageViewer({
  imagePath,
  snapshotId,
  projectId,
  datasetId,
  onZoomChange,
  onNavigateToSnapshot,
  onStructureChanged,
}: ImageViewerProps) {
  const viewerRef = useRef<HTMLDivElement>(null)
  const containerRef = useRef<HTMLDivElement>(null)

  const [zoom, setZoom] = useState(1)
  const baseZoomRef = useRef(1)
  const [position, setPosition] = useState({ x: 0, y: 0 })
  const [error, setError] = useState<string | null>(null)
  const [imageLoading, setImageLoading] = useState(false)
  const [imageReloadKey, setImageReloadKey] = useState(0)

  const [snapshotInput, setSnapshotInput] = useState<Record<string, unknown> | null>(null)
  const [snapshotOutput, setSnapshotOutput] = useState<Record<string, unknown> | null>(null)
  const [parentSnapshotId, setParentSnapshotId] = useState<string | null>(null)
  const [snapshotCreateTime, setSnapshotCreateTime] = useState<string>('')
  const [snapshotLoading, setSnapshotLoading] = useState(false)
  const [snapshotError, setSnapshotError] = useState<string | null>(null)
  const [snapshotReloadKey, setSnapshotReloadKey] = useState(0)
  const [topPanelPercent, setTopPanelPercent] = useState(50)
  const [isPanelResizing, setIsPanelResizing] = useState(false)

  const [isDragging, setIsDragging] = useState(false)
  const dragStart = useRef({ x: 0, y: 0 })
  const positionStart = useRef({ x: 0, y: 0 })

  useEffect(() => {
    setError(null)
    setImageLoading(Boolean(imagePath))
  }, [imagePath])

  const isClusteringSnapshot = Boolean(snapshotId && (/^s_cluster_/i.test(snapshotId) || /^cluster_sub_/i.test(snapshotId)))
  const isPreprocessingSnapshot = Boolean(snapshotId && /^s_pre_/i.test(snapshotId))
  const isDgeSnapshot = Boolean(snapshotId && /^s_dge_/i.test(snapshotId))
  const isAnnotationSnapshot = Boolean(snapshotId && /^annot_/i.test(snapshotId))
  const shouldShowSnapshotPanel = isClusteringSnapshot || isPreprocessingSnapshot || isDgeSnapshot || isAnnotationSnapshot

  useEffect(() => {
    let cancelled = false
    const loadSnapshotOutput = async () => {
      if (!shouldShowSnapshotPanel || !snapshotId) {
        setSnapshotInput(null)
        setSnapshotOutput(null)
        setParentSnapshotId(null)
        setSnapshotCreateTime('')
        setSnapshotError(null)
        return
      }
      setSnapshotLoading(true)
      setSnapshotError(null)
      try {
        const detail = await getSnapshotDetails(snapshotId)
        const input = detail.input && typeof detail.input === 'object'
          ? (detail.input as Record<string, unknown>)
          : {}
        const output = detail.output && typeof detail.output === 'object'
          ? (detail.output as Record<string, unknown>)
          : {}
        if (!cancelled) {
          setSnapshotInput(input)
          setSnapshotOutput(output)
          setParentSnapshotId(detail.parent_snapshot_id ?? null)
          setSnapshotCreateTime(detail.create_time || '')
        }
      } catch (e) {
        if (!cancelled) {
          setSnapshotInput(null)
          setSnapshotOutput(null)
          setParentSnapshotId(null)
          setSnapshotCreateTime('')
          setSnapshotError(e instanceof Error ? e.message : 'Failed to load snapshot details')
        }
      } finally {
        if (!cancelled) {
          setSnapshotLoading(false)
        }
      }
    }
    void loadSnapshotOutput()
    return () => {
      cancelled = true
    }
  }, [shouldShowSnapshotPanel, snapshotId, snapshotReloadKey])

  useEffect(() => {
    if (!isPanelResizing) return

    const handleMouseMove = (e: MouseEvent) => {
      const root = viewerRef.current
      if (!root) return
      const rect = root.getBoundingClientRect()
      if (rect.height <= 0) return
      const y = e.clientY - rect.top
      const percent = (y / rect.height) * 100
      const clamped = Math.min(75, Math.max(25, percent))
      setTopPanelPercent(clamped)
    }

    const handleMouseUp = () => {
      setIsPanelResizing(false)
    }

    window.addEventListener('mousemove', handleMouseMove)
    window.addEventListener('mouseup', handleMouseUp)
    return () => {
      window.removeEventListener('mousemove', handleMouseMove)
      window.removeEventListener('mouseup', handleMouseUp)
    }
  }, [isPanelResizing])

  const handleWheel = (e: React.WheelEvent) => {
    e.preventDefault()
    setZoom((z) => {
      const factor = e.deltaY > 0 ? 0.85 : 1.15
      const next = z * factor
      const newZoom = Math.min(Math.max(baseZoomRef.current * 0.5, next), 4)
      onZoomChange?.(newZoom)
      return newZoom
    })
  }

  const handleMouseDown = (e: React.MouseEvent) => {
    if (e.button !== 0) return
    setIsDragging(true)
    dragStart.current = { x: e.clientX, y: e.clientY }
    positionStart.current = position
  }

  useEffect(() => {
    if (!isDragging) return

    const handleWindowMouseMove = (e: MouseEvent) => {
      setPosition({
        x: positionStart.current.x + e.clientX - dragStart.current.x,
        y: positionStart.current.y + e.clientY - dragStart.current.y,
      })
    }

    const handleWindowMouseUp = () => setIsDragging(false)

    window.addEventListener('mousemove', handleWindowMouseMove)
    window.addEventListener('mouseup', handleWindowMouseUp)

    return () => {
      window.removeEventListener('mousemove', handleWindowMouseMove)
      window.removeEventListener('mouseup', handleWindowMouseUp)
    }
  }, [isDragging])

  const handleDoubleClick = () => {
    setZoom(baseZoomRef.current)
    setPosition({ x: 0, y: 0 })
    onZoomChange?.(baseZoomRef.current)
  }

  const handleImageLoad = (e: React.SyntheticEvent<HTMLImageElement>) => {
    const img = e.currentTarget
    const container = containerRef.current
    if (!container) return

    const cw = container.clientWidth
    const ch = container.clientHeight
    const iw = img.naturalWidth
    const ih = img.naturalHeight

    const padding = 40
    const fitScale = Math.min((cw - padding) / iw, (ch - padding) / ih)
    const initialScale = Math.min(1, fitScale * 1.15)

    baseZoomRef.current = initialScale
    setZoom(initialScale)
    setPosition({ x: 0, y: 0 })
    setImageLoading(false)
    onZoomChange?.(initialScale)
  }

  const handleImageError = () => {
    setImageLoading(false)
    setError('Failed to load image')
  }

  const handleStructureMutationSuccess = () => {
    // Force-refresh both image and snapshot detail payload for overwrite flows.
    setError(null)
    setImageLoading(Boolean(imagePath))
    setImageReloadKey((prev) => prev + 1)
    setSnapshotReloadKey((prev) => prev + 1)
    onStructureChanged?.()
  }

  if (!imagePath) {
    return (
      <div className="pdf-viewer-empty">
        <Empty
          image={<FileImageOutlined style={{ fontSize: 64, color: '#ccc' }} />}
          description="Select an image to view"
        />
      </div>
    )
  }

  const normalizedPath = imagePath.replace(/\\/g, '/')
  const posixAbsolute = normalizedPath.startsWith('/') ? normalizedPath : `/${normalizedPath}`
  const imageUrl = `local-file://${encodeURI(posixAbsolute)}?v=${imageReloadKey}`
  const summaryText = isPreprocessingSnapshot
    ? buildPreprocessingSummary(snapshotCreateTime || 'N/A', snapshotInput || {}, snapshotOutput || {})
    : buildClusteringSummary(snapshotCreateTime || 'N/A', snapshotInput || {}, snapshotOutput || {})
  const clustersText = snapshotOutput?.clusters_summary
    ? JSON.stringify(snapshotOutput.clusters_summary, null, 2)
    : 'No cluster summary available.'
  const clusterIds = (snapshotOutput?.cluster_ids && typeof snapshotOutput.cluster_ids === 'object'
    ? (snapshotOutput.cluster_ids as Record<string, string>)
    : {}) as Record<string, string>
  const clusterSummary = (snapshotOutput?.summary && typeof snapshotOutput.summary === 'object'
    ? (snapshotOutput.summary as Record<string, unknown>)
    : {}) as Record<string, unknown>
  const clusterMethod = (clusterSummary.method ?? snapshotInput?.method ?? 'leiden') as 'leiden' | 'louvain' | 'cplearn'
  const clusterResolution = typeof clusterSummary.resolution === 'number'
    ? clusterSummary.resolution
    : (typeof snapshotInput?.resolution === 'number' ? snapshotInput.resolution : undefined)
  const lineageTab = {
    key: 'lineage',
    label: 'Lineage',
    children: snapshotId ? (
      <LineagePanel
        snapshotId={snapshotId}
        onNavigateToSnapshot={onNavigateToSnapshot}
      />
    ) : null,
  }

  return (
    <div
      className={`pdf-viewer-container ${isPanelResizing ? 'panel-resizing' : ''}`}
      ref={viewerRef}
    >
      <div
        className="pdf-image-area"
        ref={containerRef}
        onWheel={handleWheel}
        onMouseDown={handleMouseDown}
        onDoubleClick={handleDoubleClick}
        style={{
          cursor: isDragging ? 'grabbing' : 'grab',
          flex: `0 0 ${topPanelPercent}%`,
          height: `${topPanelPercent}%`,
        }}
      >
        {error ? (
          <Text type="secondary">{error}</Text>
        ) : (
          <>
            {imageLoading && (
              <div className="image-loading-overlay">
                <Spin size="small" />
              </div>
            )}
            <img
              key={`${imagePath ?? 'no-image'}-${imageReloadKey}`}
              src={imageUrl}
              className="pdf-image"
              draggable={false}
              alt="view"
              onLoad={handleImageLoad}
              onError={handleImageError}
              style={{
                transform: `translate(${position.x}px, ${position.y}px) scale(${zoom})`,
              }}
            />
          </>
        )}
      </div>
      <div
        className="pdf-panel-resizer"
        onMouseDown={() => setIsPanelResizing(true)}
      />
      <div
        className="pdf-bottom-area"
        style={{
          flex: `0 0 ${100 - topPanelPercent}%`,
          height: `${100 - topPanelPercent}%`,
        }}
      >
        {shouldShowSnapshotPanel ? (
          snapshotLoading ? (
            <div className="snapshot-panel-loading">
              <Spin size="small" />
            </div>
          ) : snapshotError ? (
            <Text type="secondary">{snapshotError}</Text>
          ) : (
            <Tabs
              size="small"
              items={
                isClusteringSnapshot
                  ? [
                      {
                        key: 'summary',
                        label: 'Summary',
                        children: <div className="snapshot-tab-text">{summaryText}</div>,
                      },
                      {
                        key: 'clusters',
                        label: 'Clusters',
                        children: (
                          <div>
                            {Object.keys(clusterIds).length > 0 ? (
                              <ClusterActionsTable
                                projectId={projectId}
                                datasetId={datasetId}
                                snapshotId={snapshotId}
                                method={clusterMethod}
                                resolution={clusterResolution}
                                clusterIds={clusterIds}
                                onRunSuccess={handleStructureMutationSuccess}
                              />
                            ) : (
                              <pre className="snapshot-tab-pre">{clustersText}</pre>
                            )}
                          </div>
                        ),
                      },
                      lineageTab,
                    ]
                  : isAnnotationSnapshot
                    ? [
                        {
                          key: 'info',
                          label: 'Info',
                          children: (
                            <AnnotationInfoPanel
                              imagePath={imagePath}
                              snapshotId={snapshotId || ''}
                              parentSnapshotId={parentSnapshotId}
                              snapshotOutput={snapshotOutput}
                              onNavigateToSnapshot={onNavigateToSnapshot}
                              onRenameSuccess={() => {
                                handleStructureMutationSuccess()
                              }}
                            />
                          ),
                        },
                        lineageTab,
                      ]
                  : [
                      {
                        key: 'summary',
                        label: 'Summary',
                        children: <div className="snapshot-tab-text">{summaryText}</div>,
                      },
                      lineageTab,
                    ]
              }
            />
          )
        ) : null}
      </div>
    </div>
  )
}
