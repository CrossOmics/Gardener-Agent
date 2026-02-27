import { useState, useCallback, useEffect, useMemo, useRef } from 'react'
import { ConfigProvider, Layout, Button, Empty, Typography, message } from 'antd'
import { MenuUnfoldOutlined, MessageOutlined, LeftOutlined, RightOutlined, FolderOutlined, FolderOpenOutlined, FileImageOutlined, ProjectOutlined, DatabaseOutlined, CopyOutlined, DownloadOutlined } from '@ant-design/icons'
import ProjectSidebar from './components/ProjectSidebar'
import type { NavItem } from './components/ProjectSidebar'
import ChatPanel from './components/ChatPanel'
import AnalysisSettings from './components/AnalysisSettings'
import AnnotationTools from './components/AnnotationTools'
import ImageViewer from './components/ImageViewer'
import FolderSummary from './components/FolderSummary'
import ProjectSummary from './components/ProjectSummary'
import PipelineSetup from './components/PipelineSetup'
import DatasetStageRunner from './components/DatasetStageRunner'
import { useResizablePanels } from './hooks/useResizablePanels'
import { useArrowNavigation } from './hooks/useArrowNavigation'
import { useTreeRefreshCoordinator, type SidebarRefreshFn } from './hooks/useTreeRefreshCoordinator'
import './App.css'

const { Content } = Layout
const { Title, Text } = Typography

type ViewMode = 'chat' | 'settings' | 'annotation' | 'pipeline-setup'

interface PipelineSetupInfo {
  datasetId: string
  datasetName: string
}

function useProjectId(): string | null {
  return useMemo(() => {
    const params = new URLSearchParams(window.location.search)
    return params.get('projectId')
  }, [])
}

function useIsWelcome(): boolean {
  return useMemo(() => {
    const params = new URLSearchParams(window.location.search)
    return params.get('welcome') === 'true'
  }, [])
}

export default function App() {
  const projectId = useProjectId()
  const isWelcome = useIsWelcome()

  const [sidebarCollapsed, setSidebarCollapsed] = useState(false)
  const [chatCollapsed, setChatCollapsed] = useState(false)
  const [currentNavItem, setCurrentNavItem] = useState<NavItem | null>(null)
  const [viewMode, setViewMode] = useState<ViewMode>('chat')
  const [navList, setNavList] = useState<NavItem[]>([])
  const [imageZoom, setImageZoom] = useState<number | null>(null)
  const [pipelineSetupInfo, setPipelineSetupInfo] = useState<PipelineSetupInfo | null>(null)
  const [copyingImage, setCopyingImage] = useState(false)
  const [savingImage, setSavingImage] = useState(false)
  const [messageApi, contextHolder] = message.useMessage()

  const sidebarRefreshRef = useRef<SidebarRefreshFn | null>(null)
  const { refreshTreeKeepingPosition } = useTreeRefreshCoordinator({
    currentNavItem,
    sidebarRefreshRef,
  })

  const {
    hasPrev,
    hasNext,
    handlePrev,
    handleNext,
  } = useArrowNavigation(navList, currentNavItem, setCurrentNavItem)

  const {
    chatWidth,
    sidebarWidth,
    isResizing,
    isSidebarResizing,
    handleChatMouseDown,
    handleSidebarMouseDown,
  } = useResizablePanels()

  const handleNavigate = useCallback((item: NavItem) => {
    setCurrentNavItem((prev) => {
      if (prev && prev.type === item.type && prev.path === item.path && prev.id === item.id) {
        return prev
      }
      return item
    })
    if (viewMode === 'pipeline-setup') {
      setViewMode('chat')
      setPipelineSetupInfo(null)
    }
  }, [viewMode])

  const handleStartPipelineSetup = useCallback((datasetId: string, datasetName: string) => {
    setPipelineSetupInfo({ datasetId, datasetName })
    setViewMode('pipeline-setup')
    setCurrentNavItem(null)
  }, [])

  const handlePipelineSetupBack = useCallback(() => {
    setViewMode('chat')
    setPipelineSetupInfo(null)
  }, [])

  const handlePipelineComplete = useCallback((datasetId: string, datasetName: string) => {
    const datasetNavItem: NavItem = {
      type: 'dataset',
      id: datasetId,
      path: datasetId,
      name: datasetName,
      projectId: projectId ?? undefined,
      datasetId,
    }
    setCurrentNavItem(datasetNavItem)
    setViewMode('chat')
    setPipelineSetupInfo(null)
    void sidebarRefreshRef.current?.({ preferredNavItem: datasetNavItem })
  }, [projectId])

  const handleStageRunComplete = useCallback(() => {
    void refreshTreeKeepingPosition()
  }, [refreshTreeKeepingPosition])

  const handleNavigateToSnapshot = useCallback((targetSnapshotId: string) => {
    const folderTarget = navList.find(
      (item) => item.type === 'snapshot' && item.datasetId === currentNavItem?.datasetId && item.snapshotId === targetSnapshotId
    )
    if (folderTarget) {
      setCurrentNavItem(folderTarget)
      return
    }

    const imageTarget = navList.find(
      (item) => item.type === 'image' && item.datasetId === currentNavItem?.datasetId && item.snapshotId === targetSnapshotId
    )
    if (imageTarget) {
      setCurrentNavItem(imageTarget)
    }
  }, [currentNavItem?.datasetId, navList])

  useEffect(() => {
    setImageZoom(null)
  }, [currentNavItem?.path])

  const getParentFolder = (path: string): string | null => {
    const separator = path.includes('/') ? '/' : '\\'
    const parts = path.split(separator)
    if (parts.length >= 2) {
      return parts[parts.length - 2]
    }
    return null
  }

  const activeDatasetId = currentNavItem?.datasetId ?? null
  const activeSnapshotId = currentNavItem?.snapshotId ?? null
  const isImageNav = currentNavItem?.type === 'image'
  const currentImagePath = isImageNav ? currentNavItem.path : ''

  const handleCopyImage = useCallback(async () => {
    if (!currentImagePath) return
    if (!window.electronAPI?.copyImageToClipboard) {
      messageApi.error('Copy Image is unavailable. Rebuild Electron and restart the app.')
      return
    }
    setCopyingImage(true)
    try {
      await window.electronAPI.copyImageToClipboard(currentImagePath)
      messageApi.success('Image copied')
    } catch (error) {
      messageApi.error(error instanceof Error ? error.message : 'Failed to copy image')
    } finally {
      setCopyingImage(false)
    }
  }, [currentImagePath, messageApi])

  const handleSaveImageAs = useCallback(async () => {
    if (!currentImagePath) return
    if (!window.electronAPI?.saveImageAs) {
      messageApi.error('Save As is unavailable. Rebuild Electron and restart the app.')
      return
    }
    setSavingImage(true)
    try {
      const result = await window.electronAPI.saveImageAs(currentImagePath)
      if (!result.canceled) {
        messageApi.success('Image saved')
      }
    } catch (error) {
      messageApi.error(error instanceof Error ? error.message : 'Failed to save image')
    } finally {
      setSavingImage(false)
    }
  }, [currentImagePath, messageApi])

  if (isWelcome || !projectId) {
    return (
      <ConfigProvider
        theme={{
          token: {
            colorPrimary: '#000000',
            colorBgContainer: '#ffffff',
            colorBgElevated: '#ffffff',
            colorBgLayout: '#f5f5f5',
            colorBorder: '#e8e8e8',
            colorText: 'rgba(0, 0, 0, 0.88)',
            colorTextSecondary: 'rgba(0, 0, 0, 0.65)',
            colorTextTertiary: 'rgba(0, 0, 0, 0.45)',
            borderRadius: 8,
          },
        }}
      >
        <div className="welcome-page">
          <div className="welcome-content">
            <Title level={2}>Welcome to Lotus</Title>
            <Text type="secondary" style={{ fontSize: 16, marginBottom: 24, display: 'block' }}>
              Single-cell analysis platform
            </Text>
            <Empty
              description={
                <span>
                  No project selected.<br />
                  Use <strong>File -&gt; Open Project</strong> to open an existing project,<br />
                  or <strong>File -&gt; New Project</strong> to create a new one.
                </span>
              }
            />
          </div>
        </div>
      </ConfigProvider>
    )
  }

  return (
    <ConfigProvider
      theme={{
        token: {
          colorPrimary: '#000000',
          colorBgContainer: '#ffffff',
          colorBgElevated: '#ffffff',
          colorBgLayout: '#f5f5f5',
          colorBorder: '#e8e8e8',
          colorText: 'rgba(0, 0, 0, 0.88)',
          colorTextSecondary: 'rgba(0, 0, 0, 0.65)',
          colorTextTertiary: 'rgba(0, 0, 0, 0.45)',
          borderRadius: 8,
        },
        components: {
          Button: {
            borderRadius: 20,
            defaultBg: 'transparent',
            defaultHoverBg: '#ececec',
            defaultActiveBg: '#e0e0e0',
            defaultHoverBorderColor: 'transparent',
            defaultBorderColor: '#d9d9d9',
            primaryColor: '#ffffff',
          },
        },
      }}
    >
      {contextHolder}
      <Layout className={`layout ${isSidebarResizing || isResizing ? 'is-resizing' : ''}`}>
        {!sidebarCollapsed && (
          <>
            <div className="sidebar" style={{ width: sidebarWidth, flex: `0 0 ${sidebarWidth}px` }}>
              <ProjectSidebar
                projectId={projectId}
                onCollapse={() => setSidebarCollapsed(true)}
                currentNavItem={currentNavItem}
                onNavigate={handleNavigate}
                suspendAutoNavigate={viewMode !== 'chat'}
                onOpenSettings={() => {
                  setCurrentNavItem(null)
                  setViewMode('settings')
                }}
                onOpenAnnotation={() => {
                  setCurrentNavItem(null)
                  setViewMode('annotation')
                }}
                onNavListChange={setNavList}
                onStartPipelineSetup={handleStartPipelineSetup}
                onRefreshRef={sidebarRefreshRef}
              />
            </div>
            <div
              className={`sidebar-resizer ${isSidebarResizing ? 'resizing' : ''}`}
              onMouseDown={handleSidebarMouseDown}
            />
          </>
        )}

        {sidebarCollapsed && (
          <div className="collapsed-sidebar">
            <Button
              type="text"
              icon={<MenuUnfoldOutlined />}
              onClick={() => setSidebarCollapsed(false)}
              title="Show Sidebar"
            />
          </div>
        )}

        <Content className="content content-main">
          {currentNavItem ? (
            <div className="content-viewer">
              <div className="content-header">
                <Button
                  type="text"
                  icon={<LeftOutlined />}
                  className="nav-header-btn"
                  onClick={handlePrev}
                  disabled={!hasPrev}
                />
                <div className="content-header-title">
                  {currentNavItem.type === 'project' ? (
                    <>
                      <ProjectOutlined className="content-header-icon project" />
                      <span className="content-header-name">{currentNavItem.name}</span>
                    </>
                  ) : currentNavItem.type === 'dataset' ? (
                    <>
                      <DatabaseOutlined className="content-header-icon dataset" />
                      <span className="content-header-name">{currentNavItem.name}</span>
                    </>
                  ) : currentNavItem.type === 'image' ? (
                    <>
                      {getParentFolder(currentNavItem.path) && (
                        <>
                          <FolderOutlined className="content-header-icon folder" />
                          <span className="content-header-breadcrumb">{getParentFolder(currentNavItem.path)}</span>
                          <span className="content-header-separator">/</span>
                        </>
                      )}
                      <FileImageOutlined className="content-header-icon image" />
                      <span className="content-header-name">{currentNavItem.name}</span>
                      {imageZoom !== null && (
                        <span className="content-header-zoom">{Math.round(imageZoom * 100)}%</span>
                      )}
                    </>
                  ) : (
                    <>
                      <FolderOpenOutlined className="content-header-icon folder" />
                      <span className="content-header-name">{currentNavItem.name}</span>
                    </>
                  )}
                </div>
                {isImageNav && (
                  <div className="image-action-group">
                    <Button
                      size="small"
                      className="image-action-btn"
                      icon={<CopyOutlined />}
                      onClick={() => void handleCopyImage()}
                      loading={copyingImage}
                    >
                      Copy Image
                    </Button>
                    <Button
                      size="small"
                      className="image-action-btn"
                      icon={<DownloadOutlined />}
                      onClick={() => void handleSaveImageAs()}
                      loading={savingImage}
                    >
                      Save As
                    </Button>
                  </div>
                )}
                <Button
                  type="text"
                  icon={<RightOutlined />}
                  className="nav-header-btn"
                  onClick={handleNext}
                  disabled={!hasNext}
                />
              </div>
              {currentNavItem.type === 'project' ? (
                <ProjectSummary projectId={currentNavItem.projectId || currentNavItem.id} projectName={currentNavItem.name} />
              ) : currentNavItem.type === 'dataset' ? (
                <DatasetStageRunner
                  projectId={currentNavItem.projectId || projectId}
                  datasetId={currentNavItem.datasetId || ''}
                  datasetName={currentNavItem.name}
                  onRunSuccess={handleStageRunComplete}
                />
              ) : currentNavItem.type === 'image' ? (
                <ImageViewer
                  imagePath={currentNavItem.path}
                  snapshotId={currentNavItem.snapshotId}
                  projectId={currentNavItem.projectId}
                  datasetId={currentNavItem.datasetId}
                  onZoomChange={setImageZoom}
                  onNavigateToSnapshot={handleNavigateToSnapshot}
                  onStructureChanged={handleStageRunComplete}
                />
              ) : currentNavItem.type === 'snapshot' ? (
                <FolderSummary
                  folderPath={currentNavItem.path}
                  folderName={currentNavItem.name}
                  projectId={currentNavItem.projectId}
                  datasetId={currentNavItem.datasetId}
                  snapshotId={currentNavItem.snapshotId}
                  isSnapshotFolder={true}
                  snapshotCategory={currentNavItem.snapshotCategory}
                  onRunSuccess={handleStageRunComplete}
                />
              ) : (
                <FolderSummary
                  folderPath={currentNavItem.path}
                  folderName={currentNavItem.name}
                  projectId={currentNavItem.projectId}
                  datasetId={currentNavItem.datasetId}
                  snapshotId={currentNavItem.snapshotId}
                  isSnapshotFolder={currentNavItem.isSnapshotFolder}
                  snapshotCategory={currentNavItem.snapshotCategory}
                  onRunSuccess={handleStageRunComplete}
                />
              )}
            </div>
          ) : viewMode === 'pipeline-setup' && pipelineSetupInfo ? (
            <PipelineSetup
              projectId={projectId}
              datasetId={pipelineSetupInfo.datasetId}
              datasetName={pipelineSetupInfo.datasetName}
              onBack={handlePipelineSetupBack}
              onComplete={handlePipelineComplete}
            />
          ) : viewMode === 'settings' ? (
            <AnalysisSettings projectId={projectId} />
          ) : viewMode === 'annotation' ? (
            <AnnotationTools projectId={projectId} />
          ) : (
            <ChatPanel
              projectId={projectId}
              datasetId={activeDatasetId}
              snapshotId={activeSnapshotId}
              onAgentResponse={handleStageRunComplete}
            />
          )}
        </Content>

        {currentNavItem && !chatCollapsed && (
          <>
            <div
              className={`chat-resizer ${isResizing ? 'resizing' : ''}`}
              onMouseDown={handleChatMouseDown}
            />
            <Content className="content content-chat" style={{ width: chatWidth, flex: `0 0 ${chatWidth}px` }}>
              <ChatPanel
                projectId={projectId}
                datasetId={activeDatasetId}
                snapshotId={activeSnapshotId}
                showHeader
                onClose={() => setChatCollapsed(true)}
                onAgentResponse={handleStageRunComplete}
              />
            </Content>
          </>
        )}

        {currentNavItem && chatCollapsed && (
          <button
            className="floating-btn floating-btn-right"
            onClick={() => setChatCollapsed(false)}
            title="Show Chat"
          >
            <MessageOutlined />
          </button>
        )}
      </Layout>
    </ConfigProvider>
  )
}
