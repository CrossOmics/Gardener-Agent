import { app, BrowserWindow, Menu, dialog } from 'electron'
import { getPaths, resolveDataRoot } from './src/config'
import { createProjectService } from './src/services/projectService'
import { createWindowManager } from './src/windows/windowManager'
import { buildAppMenu } from './src/menu/buildAppMenu'
import { setupIpcHandlers } from './src/ipc/setupIpcHandlers'
import { registerLocalFileProtocol, registerLocalFileSchemeAsPrivileged } from './src/protocol/localFile'
import { createLocalServiceManager } from './src/services/localServiceManager'
import path from 'path'

const HIDDEN_FOLDERS = ['raw_data', 'snapshots', 'snapshots_anndata']

registerLocalFileSchemeAsPrivileged()

let runtime: {
  projectService: ReturnType<typeof createProjectService>
  windowManager: ReturnType<typeof createWindowManager>
  localServiceManager: ReturnType<typeof createLocalServiceManager>
} | null = null
let localServiceManagerRef: ReturnType<typeof createLocalServiceManager> | null = null
let loadingWindow: BrowserWindow | null = null
let isQuitting = false
let isStoppingServices = false

function createLoadingWindow() {
  const win = new BrowserWindow({
    width: 620,
    height: 360,
    resizable: false,
    minimizable: false,
    maximizable: false,
    frame: true,
    title: 'Gardener - Starting Services',
    show: false,
    backgroundColor: '#f6f8fb',
  })

  const html = `<!doctype html>
<html>
  <head>
    <meta charset="utf-8" />
    <title>Starting Gardener</title>
    <style>
      body { margin: 0; font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif; background: #f6f8fb; color: #111827; }
      .wrap { height: 100vh; display: flex; flex-direction: column; align-items: center; justify-content: center; gap: 14px; padding: 24px; box-sizing: border-box; }
      .title { font-size: 20px; font-weight: 600; }
      .status { font-size: 14px; color: #374151; text-align: center; max-width: 520px; line-height: 1.5; }
      .hint { font-size: 12px; color: #6b7280; }
      .spinner { width: 28px; height: 28px; border-radius: 50%; border: 3px solid #d1d5db; border-top-color: #111827; animation: spin 0.9s linear infinite; }
      @keyframes spin { to { transform: rotate(360deg); } }
    </style>
  </head>
  <body>
    <div class="wrap">
      <div class="spinner"></div>
      <div class="title">Starting Gardener</div>
      <div class="status" id="status">Initializing services...</div>
      <div class="hint">The app will open when backend is ready.</div>
    </div>
    <script>
      window.__setStatus = function (text) {
        const el = document.getElementById('status');
        if (el) el.textContent = text;
      };
    </script>
  </body>
</html>`

  void win.loadURL(`data:text/html;charset=UTF-8,${encodeURIComponent(html)}`)
  win.once('ready-to-show', () => win.show())
  return win
}

function updateLoadingStatus(text: string) {
  if (!loadingWindow || loadingWindow.isDestroyed()) return
  const safe = JSON.stringify(text)
  void loadingWindow.webContents.executeJavaScript(`window.__setStatus && window.__setStatus(${safe})`)
}

app.whenReady().then(() => {
  const bootstrap = async () => {
    loadingWindow = createLoadingWindow()

    const localServiceManager = createLocalServiceManager()
    localServiceManagerRef = localServiceManager
    updateLoadingStatus('Starting backend and agent services...')
    const agentStartupPromise = localServiceManager.startService('agent')
      .then(() => {
        console.log('[local-service] agent startup completed')
      })
      .catch((error) => {
        const msg = error instanceof Error ? error.message : 'Unknown error'
        console.error(`[local-service] agent startup failed: ${msg}`)
      })

    await localServiceManager.startService('backend')
    updateLoadingStatus('Backend is ready. Loading workspace...')

    const dataRoot = await resolveDataRoot()
    const { PROJECTS_PATH, PRELOAD_PATH } = getPaths(dataRoot, __dirname)
    const rendererDevUrl = process.env.LOTUS_RENDERER_URL || 'http://localhost:5173'
    const rendererProdIndexPath = path.join(process.resourcesPath, 'react', 'dist', 'index.html')
    const projectService = createProjectService(PROJECTS_PATH, HIDDEN_FOLDERS)
    const windowManager = createWindowManager({
      preloadPath: PRELOAD_PATH,
      getProjects: projectService.getProjects,
      rendererDevUrl,
      rendererProdIndexPath,
      isPackaged: app.isPackaged,
    })
    runtime = { projectService, windowManager, localServiceManager }

    const rebuildMenu = () => {
      Menu.setApplicationMenu(buildAppMenu({
        getProjects: projectService.getProjects,
        createProjectViaApi: projectService.createProjectViaApi,
        openProjectFromCurrentContext: windowManager.openProjectFromCurrentContext,
        createProjectWindow: windowManager.createProjectWindow,
      }))
    }

    registerLocalFileProtocol()

    setupIpcHandlers({
      getProjects: projectService.getProjects,
      getProjectInfo: projectService.getProjectInfo,
      getDatasets: projectService.getDatasets,
      getDatasetFiles: projectService.getDatasetFiles,
      createProjectViaApi: projectService.createProjectViaApi,
      getWorkspacePath: projectService.getWorkspacePath,
      createProjectWindow: windowManager.createProjectWindow,
      rebuildMenu,
    })

    rebuildMenu()

    void agentStartupPromise

    const projects = projectService.getProjects()
    if (projects.length > 0) {
      windowManager.createProjectWindow(projects[0].id)
    } else {
      windowManager.createWelcomeWindow()
    }
    if (loadingWindow && !loadingWindow.isDestroyed()) {
      loadingWindow.close()
    }
    loadingWindow = null
  }

  void bootstrap().catch(async (error) => {
    const message = error instanceof Error ? error.message : 'Failed to resolve project root path'
    dialog.showErrorBox('Startup Failed', message)
    if (loadingWindow && !loadingWindow.isDestroyed()) {
      loadingWindow.close()
    }
    loadingWindow = null
    if (localServiceManagerRef) await localServiceManagerRef.stopAll()
    app.quit()
  })
})

app.on('before-quit', (event) => {
  if (isQuitting || isStoppingServices) return
  event.preventDefault()
  isStoppingServices = true
  void (async () => {
    try {
      if (localServiceManagerRef) {
        await localServiceManagerRef.stopAll()
      }
    } finally {
      isQuitting = true
      isStoppingServices = false
      app.quit()
    }
  })()
})

app.on('window-all-closed', () => {
  if (process.platform !== 'darwin') {
    app.quit()
  }
})

app.on('activate', () => {
  if (BrowserWindow.getAllWindows().length === 0 && runtime) {
    const projects = runtime.projectService.getProjects()
    if (projects.length > 0) {
      runtime.windowManager.createProjectWindow(projects[0].id)
    } else {
      runtime.windowManager.createWelcomeWindow()
    }
  }
})
