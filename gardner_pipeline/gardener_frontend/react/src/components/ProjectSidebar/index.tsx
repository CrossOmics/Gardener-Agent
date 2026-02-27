import { useState, useEffect, useCallback, useRef, type MutableRefObject } from 'react'
import { Button, Typography, Spin, message, Modal, Input } from 'antd'
import {
  UploadOutlined,
  MenuFoldOutlined,
  FolderOutlined,
  FolderOpenOutlined,
  FileImageOutlined,
  BulbOutlined,
  HighlightOutlined,
  ReloadOutlined,
  DownOutlined,
  RightOutlined,
  DatabaseOutlined,
} from '@ant-design/icons'
import { getProjectDisplayMap, importDataset } from '../../services/endpoints'
import {
  deleteDataset,
  deleteSnapshot,
  deleteSnapshotsByStage,
  renameEntityDisplayName,
} from '../../services/endpoints'
import {
  collectNavItems,
  inferSnapshotIdFromDirectoryPath,
  inferSnapshotIdFromName,
  inferSnapshotIdFromPath,
  sortFiles,
  type FileNode,
  type NavItem,
} from './fileTree'
import TreeItemActions from './TreeItemActions'
import './styles.css'

const { Text } = Typography

interface DatasetInfo {
  id: string
  name: string
  displayName?: string
  path: string
  createdAt: Date
}

interface ProjectInfo {
  id: string
  name: string
  path: string
  createdAt: Date
}

interface DatasetWithFiles extends DatasetInfo {
  files: FileNode[]
}
export type { NavItem } from './fileTree'

type StageName = 'Preprocessing' | 'Clustering' | 'DGE' | 'Annotation'

const STAGE_DISPLAY_TO_API: Record<string, StageName> = {
  preprocessing: 'Preprocessing',
  clustering: 'Clustering',
  dge: 'DGE',
  annotation: 'Annotation',
}

function inferSnapshotCategoryFromPath(path: string): NavItem['snapshotCategory'] | undefined {
  const normalized = path.replace(/\\/g, '/').toLowerCase()
  if (normalized.includes('/preprocessing/')) return 'preprocessing'
  if (normalized.includes('/clustering/')) return 'clustering'
  if (normalized.includes('/dge/')) return 'dge'
  if (normalized.includes('/annotation/')) return 'annotation'
  return undefined
}

function hasDirectImageChild(node: FileNode): boolean {
  return (node.children || []).some((child) => !child.isDirectory && child.name.toLowerCase().endsWith('.png'))
}

function findSnapshotDisplayNameInDataset(
  snapshotsByStage: Record<string, Record<string, string>> | undefined,
  snapshotId: string
): string | undefined {
  if (!snapshotsByStage) return undefined
  for (const stageMap of Object.values(snapshotsByStage)) {
    const displayName = stageMap[snapshotId]
    if (displayName) return displayName
  }
  return undefined
}

function applyDisplayNamesToFileTree(
  files: FileNode[],
  snapshotsByStage: Record<string, Record<string, string>> | undefined
): FileNode[] {
  return files.map((file) => {
    if (!file.isDirectory) return file
    const nextChildren = file.children ? applyDisplayNamesToFileTree(file.children, snapshotsByStage) : undefined
    const isSnapshotFolder = hasDirectImageChild(file)
    const snapshotId = isSnapshotFolder
      ? (inferSnapshotIdFromDirectoryPath(file.path) ?? inferSnapshotIdFromName(file.name))
      : undefined
    const displayName = snapshotId ? findSnapshotDisplayNameInDataset(snapshotsByStage, snapshotId) : undefined
    return {
      ...file,
      displayName: displayName ?? file.displayName,
      children: nextChildren,
    }
  })
}

function detectStageNameFromFolder(folder: FileNode): StageName | null {
  const normalized = folder.name.trim().toLowerCase()
  return STAGE_DISPLAY_TO_API[normalized] ?? null
}

function updateDatasetDisplayName(files: DatasetWithFiles[], datasetId: string, newName: string): DatasetWithFiles[] {
  return files.map((dataset) => (dataset.id === datasetId ? { ...dataset, displayName: newName } : dataset))
}

function updateSnapshotDisplayNameInTree(files: FileNode[], snapshotId: string, newName: string): FileNode[] {
  return files.map((file) => {
    if (!file.isDirectory) return file
    const currentSnapshotId = inferSnapshotIdFromDirectoryPath(file.path) ?? inferSnapshotIdFromName(file.name)
    const nextChildren = file.children ? updateSnapshotDisplayNameInTree(file.children, snapshotId, newName) : undefined
    if (currentSnapshotId !== snapshotId) {
      return nextChildren ? { ...file, children: nextChildren } : file
    }
    return {
      ...file,
      displayName: newName,
      children: nextChildren,
    }
  })
}

interface ProjectSidebarProps {
  projectId: string
  onCollapse: () => void
  currentNavItem: NavItem | null
  onNavigate: (item: NavItem) => void
  onOpenSettings: () => void
  onOpenAnnotation: () => void
  onNavListChange?: (navList: NavItem[]) => void
  onStartPipelineSetup?: (datasetId: string, datasetName: string) => void
  onRefreshRef?: MutableRefObject<((options?: { preferredNavItem?: NavItem | null }) => Promise<void>) | null>
  suspendAutoNavigate?: boolean
}

export default function ProjectSidebar({
  projectId,
  onCollapse,
  currentNavItem,
  onNavigate,
  onOpenSettings,
  onOpenAnnotation,
  onNavListChange,
  onStartPipelineSetup,
  onRefreshRef,
  suspendAutoNavigate = false,
}: ProjectSidebarProps) {
  const [projectInfo, setProjectInfo] = useState<ProjectInfo | null>(null)
  const [datasets, setDatasets] = useState<DatasetWithFiles[]>([])
  const [loading, setLoading] = useState(true)
  const [expandedDatasets, setExpandedDatasets] = useState<Set<string>>(new Set())
  const [expandedFolders, setExpandedFolders] = useState<Set<string>>(new Set())
  const [uploading, setUploading] = useState(false)
  const [importDisabledReason, setImportDisabledReason] = useState<string | null>(null)
  const [renameOpen, setRenameOpen] = useState(false)
  const [renameValue, setRenameValue] = useState('')
  const [renameSubmitting, setRenameSubmitting] = useState(false)
  const [renameTarget, setRenameTarget] = useState<
    | { type: 'dataset'; datasetId: string }
    | { type: 'snapshot'; datasetId: string; snapshotId: string }
    | null
  >(null)
  const currentNavItemRef = useRef<NavItem | null>(currentNavItem)

  const isSameNavItem = useCallback((a: NavItem | null, b: NavItem | null) => {
    if (!a || !b) return false
    return a.type === b.type && a.path === b.path && a.id === b.id
  }, [])

  useEffect(() => {
    currentNavItemRef.current = currentNavItem
  }, [currentNavItem])

  const findBestNavMatch = useCallback((items: NavItem[], targetNavItem: NavItem | null): NavItem | null => {
    if (!targetNavItem) return null

    const exact = items.find((item) => item.type === targetNavItem.type && item.path === targetNavItem.path)
    if (exact) return exact

    const byId = items.find(
      (item) =>
        item.type === targetNavItem.type &&
        item.id === targetNavItem.id &&
        item.datasetId === targetNavItem.datasetId
    )
    if (byId) return byId

    // Fallback to closest existing parent folder in same dataset.
    const sameDatasetFolders = items.filter(
      (item) => (item.type === 'folder' || item.type === 'snapshot') && item.datasetId === targetNavItem.datasetId
    )
    const normalizedCurrent = targetNavItem.path.replace(/\\/g, '/')
    const parentCandidate = sameDatasetFolders
      .filter((item) => normalizedCurrent.startsWith(item.path.replace(/\\/g, '/')))
      .sort((a, b) => b.path.length - a.path.length)[0]

    return parentCandidate || null
  }, [])

  const loadProjectData = useCallback(async (options?: { preferredNavItem?: NavItem | null }) => {
    if (!projectId) return

    setLoading(true)
    try {
      if (window.electronAPI) {
        // Get project info
        const info = await window.electronAPI.getProjectInfo(projectId)
        setProjectInfo(info)

        // Get datasets for this project
        const datasetList = await window.electronAPI.getDatasets(projectId)
        let projectDisplayMap: Awaited<ReturnType<typeof getProjectDisplayMap>> | null = null
        try {
          projectDisplayMap = await getProjectDisplayMap(projectId)
        } catch (error) {
          console.warn('Failed to load project display map, fallback to IDs:', error)
        }

        // Load file tree for each dataset
        const datasetsWithFiles: DatasetWithFiles[] = await Promise.all(
          datasetList.map(async (dataset: DatasetInfo) => {
            const files = await window.electronAPI.getDatasetFiles(projectId, dataset.id)
            const datasetDisplay = projectDisplayMap?.datasets?.[dataset.id]
            return {
              ...dataset,
              displayName: datasetDisplay?.dataset_name ?? dataset.displayName,
              files: applyDisplayNamesToFileTree(files, datasetDisplay?.snapshots),
            }
          })
        )

        // Keep all datasets visible even if they currently have no PNG outputs.
        setDatasets(datasetsWithFiles)

        // Collect all navigable items. Include dataset entries so navigation
        // can reliably preserve/restore dataset-level selection.
        const allNavItems = datasetsWithFiles.flatMap((dataset) => ([
          {
            type: 'dataset' as const,
            id: dataset.id,
            path: dataset.path,
            name: dataset.displayName ?? dataset.name ?? dataset.id,
            projectId,
            datasetId: dataset.id,
          },
          ...collectNavItems(dataset.files, projectId, dataset.id),
        ]))
        onNavListChange?.(allNavItems)

        const currentNav = options?.preferredNavItem ?? currentNavItemRef.current
        const preferred = findBestNavMatch(allNavItems, currentNav)
        if (preferred) {
          if (preferred.datasetId) {
            setExpandedDatasets((prev) => {
              const next = new Set(prev)
              next.add(preferred.datasetId!)
              return next
            })
          }
          if (preferred.type === 'folder' || preferred.type === 'snapshot') {
            setExpandedFolders((prev) => {
              const next = new Set(prev)
              next.add(preferred.path)
              return next
            })
          }
          if (!suspendAutoNavigate && !isSameNavItem(currentNav, preferred)) {
            onNavigate(preferred)
          }
          return
        }

        // First-load fallback: auto-expand first dataset and select first item.
        if (datasetsWithFiles.length > 0) {
          const firstDataset = datasetsWithFiles[0]
          setExpandedDatasets((prev) => {
            const next = new Set(prev)
            next.add(firstDataset.id)
            return next
          })

          if (!suspendAutoNavigate && !currentNav && allNavItems.length > 0) {
            const firstItem = allNavItems[0]
            if (firstItem.type === 'folder' || firstItem.type === 'snapshot') {
              setExpandedFolders((prev) => {
                const next = new Set(prev)
                next.add(firstItem.path)
                return next
              })
            }
            if (!isSameNavItem(currentNav, firstItem)) {
              onNavigate(firstItem)
            }
          }
        }
      }
    } catch (error) {
      console.error('Failed to load project data:', error)
      message.error('Failed to load project data')
    } finally {
      setLoading(false)
    }
  }, [projectId, onNavListChange, findBestNavMatch, onNavigate, isSameNavItem, suspendAutoNavigate])

  useEffect(() => {
    loadProjectData()
  }, [loadProjectData])

  useEffect(() => {
    setImportDisabledReason(null)
  }, [projectId])

  // Expose refresh function to parent
  useEffect(() => {
    if (onRefreshRef) {
      onRefreshRef.current = (options) => loadProjectData(options)
    }
    return () => {
      if (onRefreshRef) {
        onRefreshRef.current = null
      }
    }
  }, [onRefreshRef, loadProjectData])

  // Auto-expand parent folders when navigation changes
  useEffect(() => {
    if (!currentNavItem) return

    const path = currentNavItem.path
    const separator = path.includes('/') ? '/' : '\\'
    const parts = path.split(separator)

    const parentPaths: string[] = []
    for (let i = 1; i < parts.length; i++) {
      parentPaths.push(parts.slice(0, i).join(separator))
    }

    if (parentPaths.length > 0) {
      setExpandedFolders((prev) => {
        const newSet = new Set(prev)
        parentPaths.forEach((p) => newSet.add(p))
        return newSet
      })
    }
  }, [currentNavItem])

  const toggleDataset = (datasetId: string) => {
    setExpandedDatasets((prev) => {
      const newSet = new Set(prev)
      if (newSet.has(datasetId)) {
        newSet.delete(datasetId)
      } else {
        newSet.add(datasetId)
      }
      return newSet
    })
  }

  const handleDatasetClick = (dataset: DatasetWithFiles) => {
    toggleDataset(dataset.id)
    onNavigate({
      type: 'dataset',
      id: dataset.id,
      path: dataset.path,
      name: dataset.displayName ?? dataset.name ?? dataset.id,
      projectId,
      datasetId: dataset.id,
    })
  }

  const toggleFolder = (folder: FileNode, dsId: string) => {
    const folderPath = folder.path
    const folderName = folder.name
    const isSnapshotFolder = hasDirectImageChild(folder)
    const snapshotCategory = isSnapshotFolder ? inferSnapshotCategoryFromPath(folderPath) : undefined
    const snapshotId = isSnapshotFolder
      ? (inferSnapshotIdFromDirectoryPath(folderPath) ?? inferSnapshotIdFromName(folderName))
      : undefined
    const displayName = folder.displayName ?? folderName

    setExpandedFolders((prev) => {
      const newSet = new Set(prev)
      if (newSet.has(folderPath)) {
        newSet.delete(folderPath)
      } else {
        newSet.add(folderPath)
      }
      return newSet
    })
    onNavigate({
      type: isSnapshotFolder ? 'snapshot' : 'folder',
      id: folderName,
      path: folderPath,
      name: displayName,
      snapshotName: isSnapshotFolder ? displayName : undefined,
      projectId,
      datasetId: dsId,
      snapshotId,
      isSnapshotFolder,
      snapshotCategory,
    })
  }

  const handleImageClick = (file: FileNode, dsId: string) => {
    onNavigate({
      type: 'image',
      id: file.id,
      path: file.path,
      name: file.name,
      projectId,
      datasetId: dsId,
      snapshotId: inferSnapshotIdFromPath(file.path),
    })
  }

  const patchCurrentNavName = useCallback((matcher: (item: NavItem) => boolean, nextName: string) => {
    if (!currentNavItemRef.current) return
    if (!matcher(currentNavItemRef.current)) return
    const nextItem: NavItem = {
      ...currentNavItemRef.current,
      name: nextName,
      snapshotName: currentNavItemRef.current.type === 'snapshot' ? nextName : currentNavItemRef.current.snapshotName,
    }
    onNavigate(nextItem)
  }, [onNavigate])

  const openRenameModal = useCallback((target: { type: 'dataset'; datasetId: string } | {
    type: 'snapshot'
    datasetId: string
    snapshotId: string
  }, initialValue: string) => {
    setRenameTarget(target)
    setRenameValue(initialValue)
    setRenameOpen(true)
  }, [])

  const submitRename = useCallback(async () => {
    const trimmed = renameValue.trim()
    if (!renameTarget || !trimmed) {
      message.error('Name cannot be empty')
      return
    }

    try {
      setRenameSubmitting(true)
      if (renameTarget.type === 'dataset') {
        await renameEntityDisplayName({
          id_type: 'dataset',
          current_id: renameTarget.datasetId,
          new_name: trimmed,
        })
        setDatasets((prev) => updateDatasetDisplayName(prev, renameTarget.datasetId, trimmed))
        patchCurrentNavName(
          (item) => item.type === 'dataset' && item.datasetId === renameTarget.datasetId,
          trimmed
        )
      } else {
        await renameEntityDisplayName({
          id_type: 'snapshot',
          current_id: renameTarget.snapshotId,
          new_name: trimmed,
        })
        setDatasets((prev) =>
          prev.map((dataset) => {
            if (dataset.id !== renameTarget.datasetId) return dataset
            return {
              ...dataset,
              files: updateSnapshotDisplayNameInTree(dataset.files, renameTarget.snapshotId, trimmed),
            }
          })
        )
        patchCurrentNavName(
          (item) => item.snapshotId === renameTarget.snapshotId,
          trimmed
        )
      }

      setRenameOpen(false)
      setRenameTarget(null)
      message.success('Renamed successfully')
    } catch (error) {
      message.error(error instanceof Error ? error.message : 'Rename failed')
    } finally {
      setRenameSubmitting(false)
    }
  }, [patchCurrentNavName, renameTarget, renameValue])

  const handleDatasetDelete = useCallback((dataset: DatasetWithFiles) => {
    Modal.confirm({
      title: 'Delete Dataset',
      content: `Delete dataset "${dataset.displayName ?? dataset.name ?? dataset.id}" and all related files?`,
      okButtonProps: { danger: true },
      onOk: async () => {
        await deleteDataset(dataset.id, false)
        message.success('Dataset deleted')
        await loadProjectData({ preferredNavItem: currentNavItemRef.current })
      },
    })
  }, [loadProjectData])

  const handleDatasetKeepFinal = useCallback((dataset: DatasetWithFiles) => {
    Modal.confirm({
      title: 'Keep Final Result',
      content:
        `Clean dataset "${dataset.displayName ?? dataset.name ?? dataset.id}" and keep only the latest Annotation result when available?`,
      onOk: async () => {
        await deleteDataset(dataset.id, true)
        message.success('Dataset cleaned and final result kept')
        await loadProjectData({ preferredNavItem: currentNavItemRef.current })
      },
    })
  }, [loadProjectData])

  const handleSnapshotDelete = useCallback((snapshotId: string, snapshotDisplayName: string) => {
    Modal.confirm({
      title: 'Delete Snapshot',
      content: `Delete snapshot "${snapshotDisplayName}" and all related files?`,
      okButtonProps: { danger: true },
      onOk: async () => {
        await deleteSnapshot(snapshotId)
        message.success('Snapshot deleted')
        await loadProjectData({ preferredNavItem: currentNavItemRef.current })
      },
    })
  }, [loadProjectData])

  const handleStageDelete = useCallback((datasetId: string, stageName: StageName, keepLatest: boolean) => {
    Modal.confirm({
      title: keepLatest ? 'Keep Latest Result' : 'Delete Stage Results',
      content: keepLatest
        ? `Keep only the latest snapshot in ${stageName} and delete the rest?`
        : `Delete all snapshots in ${stageName}?`,
      okButtonProps: keepLatest ? undefined : { danger: true },
      onOk: async () => {
        await deleteSnapshotsByStage({
          dataset_id: datasetId,
          stage_name: stageName,
          keep_latest: keepLatest,
        })
        message.success(keepLatest ? 'Kept latest result' : 'Stage snapshots deleted')
        await loadProjectData({ preferredNavItem: currentNavItemRef.current })
      },
    })
  }, [loadProjectData])

  const handleUploadDataset = async () => {
    if (!window.electronAPI) {
      message.error('Electron API not available')
      return
    }

    try {
      // Open file dialog using Electron
      const result = await window.electronAPI.showUploadDialog()

      if (result.canceled || result.filePaths.length === 0) {
        return
      }

      const filePath = result.filePaths[0]
      const fileName = filePath.split(/[/\\]/).pop() || 'dataset'
      const datasetName = fileName.replace(/\.(h5ad|csv|mtx)$/i, '')

      setUploading(true)
      message.loading({ content: 'Importing dataset...', key: 'upload' })

      // Call backend API to import dataset
      const response = await importDataset({
        project_id: projectId,
        local_file_path: filePath,
        dataset_name: datasetName,
      })

      message.success({ content: `Dataset "${response.dataset_name}" imported successfully`, key: 'upload' })

      // Navigate to pipeline setup page
      if (onStartPipelineSetup) {
        onStartPipelineSetup(response.dataset_id, response.dataset_name)
      }
    } catch (error) {
      console.error('Failed to import dataset:', error)
      const errorMessage = error instanceof Error ? error.message : 'Failed to import dataset'
      const lowered = errorMessage.toLowerCase()
      if (lowered.includes('does not exist') && lowered.includes('project')) {
        const disabledText = 'Current project is not registered in backend. Please create a new project from File -> New Project.'
        setImportDisabledReason(disabledText)
        message.error({ content: disabledText, key: 'upload' })
      } else {
        message.error({ content: errorMessage, key: 'upload' })
      }
    } finally {
      setUploading(false)
    }
  }

  const renderFiles = (files: FileNode[], dsId: string, level: number = 0) => {
    // Sort files before rendering
    const sortedFiles = sortFiles(files)
    return sortedFiles.map((file) => {
      if (file.isDirectory) {
        const isExpanded = expandedFolders.has(file.path)
        const isSelected = currentNavItem?.path === file.path
        const folderDisplayName = file.displayName ?? file.name
        const isSnapshotFolder = hasDirectImageChild(file)
        const snapshotId = isSnapshotFolder
          ? (inferSnapshotIdFromDirectoryPath(file.path) ?? inferSnapshotIdFromName(file.name))
          : undefined
        const stageName = !isSnapshotFolder ? detectStageNameFromFolder(file) : null
        return (
          <div key={file.path}>
            <div
              className={`folder-item ${isSelected ? 'selected' : ''}`}
              style={{ paddingLeft: `${(level + 2) * 16}px` }}
              onClick={() => toggleFolder(file, dsId)}
            >
              {isExpanded ? <DownOutlined className="expand-icon" /> : <RightOutlined className="expand-icon" />}
              {isExpanded ? <FolderOpenOutlined className="folder-icon" /> : <FolderOutlined className="folder-icon" />}
              <Text ellipsis className="folder-name">{folderDisplayName}</Text>
              {isSnapshotFolder && snapshotId && (
                <TreeItemActions
                  items={[
                    {
                      key: 'rename',
                      label: 'Rename',
                      onClick: () => openRenameModal({ type: 'snapshot', datasetId: dsId, snapshotId }, folderDisplayName),
                    },
                    {
                      key: 'delete',
                      label: 'Delete',
                      danger: true,
                      onClick: () => handleSnapshotDelete(snapshotId, folderDisplayName),
                    },
                  ]}
                />
              )}
              {!isSnapshotFolder && stageName && (
                <TreeItemActions
                  items={[
                    {
                      key: 'delete-stage',
                      label: 'Delete',
                      danger: true,
                      onClick: () => handleStageDelete(dsId, stageName, false),
                    },
                    {
                      key: 'keep-latest',
                      label: 'Keep Latest Result',
                      onClick: () => handleStageDelete(dsId, stageName, true),
                    },
                  ]}
                />
              )}
            </div>
            {isExpanded && file.children && (
              <div className="folder-children">
                {renderFiles(file.children, dsId, level + 1)}
              </div>
            )}
          </div>
        )
      } else {
        const isSelected = currentNavItem?.path === file.path
        return (
          <div
            key={file.path}
            className={`image-item ${isSelected ? 'selected' : ''}`}
            style={{ paddingLeft: `${(level + 2) * 16 + 16}px` }}
            onClick={() => handleImageClick(file, dsId)}
          >
            <FileImageOutlined className="image-icon" />
            <Text ellipsis className="image-name">{file.name}</Text>
          </div>
        )
      }
    })
  }

  return (
    <div className="sidebar-container">
      <div className="sidebar-header">
        <Button
          icon={<UploadOutlined />}
          onClick={handleUploadDataset}
          loading={uploading}
          disabled={!!importDisabledReason}
          className="new-project-btn"
          title={importDisabledReason || 'Import Dataset'}
        >
          Import Dataset
        </Button>
        <Button
          type="text"
          icon={<ReloadOutlined />}
          onClick={() => {
            void loadProjectData()
          }}
          className="refresh-btn"
          title="Refresh"
        />
        <Button type="text" icon={<MenuFoldOutlined />} onClick={onCollapse} className="collapse-btn" />
      </div>

      {/* Project title */}
      {projectInfo && (
        <div className="project-header">
          <Text strong ellipsis className="project-title">
            {projectInfo.name}
          </Text>
        </div>
      )}

      <div className="project-list">
        {loading ? (
          <div className="loading-container">
            <Spin size="small" />
          </div>
        ) : datasets.length === 0 ? (
          <div className="empty-projects">
            <Text className="empty-title">No datasets found</Text>
            <Text className="empty-subtitle" style={{ fontSize: 12, marginTop: 8 }}>
              Click "Import Dataset" to add data
            </Text>
          </div>
        ) : (
          datasets.map((dataset) => {
            const isDatasetExpanded = expandedDatasets.has(dataset.id)
            const datasetDisplayName = dataset.displayName ?? dataset.name ?? dataset.id
            return (
              <div key={dataset.id} className="dataset-tree-item">
                {/* Dataset level */}
                <div
                  className={`dataset-item ${
                    currentNavItem?.type === 'dataset' && currentNavItem?.datasetId === dataset.id ? 'selected' : ''
                  }`}
                  style={{ paddingLeft: '16px' }}
                  onClick={() => handleDatasetClick(dataset)}
                >
                  {isDatasetExpanded ? <DownOutlined className="expand-icon" /> : <RightOutlined className="expand-icon" />}
                  <DatabaseOutlined className="dataset-icon" />
                  <Text ellipsis className="dataset-name">{datasetDisplayName}</Text>
                  <TreeItemActions
                    items={[
                      {
                        key: 'rename-dataset',
                        label: 'Rename',
                        onClick: () => openRenameModal({ type: 'dataset', datasetId: dataset.id }, datasetDisplayName),
                      },
                      {
                        key: 'delete-dataset',
                        label: 'Delete',
                        danger: true,
                        onClick: () => handleDatasetDelete(dataset),
                      },
                      {
                        key: 'keep-final',
                        label: 'Keep Final Result',
                        onClick: () => handleDatasetKeepFinal(dataset),
                      },
                    ]}
                  />
                </div>

                {/* File tree */}
                {isDatasetExpanded && (
                  <div className="dataset-children">
                    {renderFiles(dataset.files, dataset.id)}
                  </div>
                )}
              </div>
            )
          })
        )}
      </div>

      <Modal
        title="Rename"
        open={renameOpen}
        onCancel={() => {
          if (renameSubmitting) return
          setRenameOpen(false)
          setRenameTarget(null)
        }}
        onOk={() => {
          void submitRename()
        }}
        okText="Save"
        confirmLoading={renameSubmitting}
      >
        <Input
          value={renameValue}
          onChange={(e) => setRenameValue(e.target.value)}
          onPressEnter={() => {
            void submitRename()
          }}
          autoFocus
          maxLength={128}
        />
      </Modal>

      <div className="sidebar-footer">
        <div className="settings-item" onClick={onOpenSettings}>
          <BulbOutlined className="settings-icon" />
          <Text className="settings-text">Customize Analysis</Text>
        </div>
        <div className="settings-item" onClick={onOpenAnnotation}>
          <HighlightOutlined className="settings-icon" />
          <Text className="settings-text">Annotation Tools</Text>
        </div>
      </div>
    </div>
  )
}
