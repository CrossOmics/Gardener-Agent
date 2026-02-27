import { clipboard, dialog, ipcMain, nativeImage } from 'electron'
import fs from 'node:fs'
import path from 'node:path'
import type { DatasetInfo, FileNode, ProjectInfo } from '../types'

interface SetupIpcHandlersDeps {
  getProjects: () => ProjectInfo[]
  getProjectInfo: (projectId: string) => ProjectInfo | null
  getDatasets: (projectId: string) => DatasetInfo[]
  getDatasetFiles: (projectId: string, datasetId: string) => FileNode[]
  createProjectViaApi: (projectName?: string) => Promise<ProjectInfo | null>
  getWorkspacePath: () => string
  createProjectWindow: (projectId: string) => void
  rebuildMenu: () => void
}

export function setupIpcHandlers(deps: SetupIpcHandlersDeps) {
  ipcMain.handle('get-projects', () => deps.getProjects())

  ipcMain.handle('get-datasets', (_event, projectId: string) => {
    return deps.getDatasets(projectId)
  })

  ipcMain.handle('get-dataset-files', (_event, projectId: string, datasetId: string) => {
    return deps.getDatasetFiles(projectId, datasetId)
  })

  ipcMain.handle('get-workspace-path', () => deps.getWorkspacePath())

  ipcMain.handle('get-project-info', (_event, projectId: string) => {
    return deps.getProjectInfo(projectId)
  })

  ipcMain.handle('create-project', async (_event, projectName?: string) => {
    const project = await deps.createProjectViaApi(projectName)
    if (project) {
      deps.rebuildMenu()
    }
    return project
  })

  ipcMain.handle('open-project-window', (_event, projectId: string) => {
    deps.createProjectWindow(projectId)
    return true
  })

  ipcMain.handle('show-upload-dialog', async () => {
    return dialog.showOpenDialog({
      title: 'Select Dataset File',
      filters: [
        { name: 'Supported Dataset Files', extensions: ['h5ad', 'csv', 'mtx'] },
        { name: 'AnnData Files', extensions: ['h5ad'] },
        { name: 'CSV Files', extensions: ['csv'] },
        { name: 'MTX Files', extensions: ['mtx'] },
        { name: 'All Files', extensions: ['*'] },
      ],
      properties: ['openFile'],
    })
  })

  ipcMain.handle('copy-image-to-clipboard', async (_event, imagePath: string) => {
    const sourcePath = String(imagePath || '').trim()
    if (!sourcePath) {
      throw new Error('Image path is required')
    }
    const image = nativeImage.createFromPath(sourcePath)
    if (image.isEmpty()) {
      throw new Error('Failed to load image for clipboard')
    }
    clipboard.writeImage(image)
    return true
  })

  ipcMain.handle('save-image-as', async (_event, imagePath: string) => {
    const sourcePath = String(imagePath || '').trim()
    if (!sourcePath) {
      throw new Error('Image path is required')
    }
    if (!fs.existsSync(sourcePath)) {
      throw new Error('Source image does not exist')
    }

    const ext = path.extname(sourcePath) || '.png'
    const defaultName = path.basename(sourcePath)
    const result = await dialog.showSaveDialog({
      title: 'Save Image As',
      defaultPath: defaultName,
      filters: [
        { name: ext.replace('.', '').toUpperCase() + ' Image', extensions: [ext.replace('.', '')] },
        { name: 'All Files', extensions: ['*'] },
      ],
    })

    if (result.canceled || !result.filePath) {
      return { canceled: true }
    }

    await fs.promises.copyFile(sourcePath, result.filePath)
    return { canceled: false, filePath: result.filePath }
  })
}
