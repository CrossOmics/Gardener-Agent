import { contextBridge, ipcRenderer } from 'electron'

export interface ProjectInfo {
  id: string
  name: string
  path: string
  createdAt: Date
}

export interface DatasetInfo {
  id: string
  name: string
  path: string
  createdAt: Date
}

export interface FileNode {
  id: string
  name: string
  path: string
  isDirectory: boolean
  children?: FileNode[]
}

export interface DialogResult {
  canceled: boolean
  filePaths: string[]
}

export interface SaveImageResult {
  canceled: boolean
  filePath?: string
}

contextBridge.exposeInMainWorld('electronAPI', {
  // Project operations
  getProjects: () => ipcRenderer.invoke('get-projects'),
  getProjectInfo: (projectId: string) => ipcRenderer.invoke('get-project-info', projectId),
  createProject: (projectName?: string) => ipcRenderer.invoke('create-project', projectName),
  openProjectWindow: (projectId: string) => ipcRenderer.invoke('open-project-window', projectId),

  // Dataset operations
  getDatasets: (projectId: string) => ipcRenderer.invoke('get-datasets', projectId),
  getDatasetFiles: (projectId: string, datasetId: string) =>
    ipcRenderer.invoke('get-dataset-files', projectId, datasetId),

  // File dialogs
  showUploadDialog: () => ipcRenderer.invoke('show-upload-dialog'),
  copyImageToClipboard: (imagePath: string) => ipcRenderer.invoke('copy-image-to-clipboard', imagePath),
  saveImageAs: (imagePath: string) => ipcRenderer.invoke('save-image-as', imagePath),

  // Workspace
  getWorkspacePath: () => ipcRenderer.invoke('get-workspace-path'),
})

export interface ElectronAPI {
  // Project operations
  getProjects: () => Promise<ProjectInfo[]>
  getProjectInfo: (projectId: string) => Promise<ProjectInfo | null>
  createProject: (projectName?: string) => Promise<ProjectInfo | null>
  openProjectWindow: (projectId: string) => Promise<boolean>

  // Dataset operations
  getDatasets: (projectId: string) => Promise<DatasetInfo[]>
  getDatasetFiles: (projectId: string, datasetId: string) => Promise<FileNode[]>

  // File dialogs
  showUploadDialog: () => Promise<DialogResult>
  copyImageToClipboard: (imagePath: string) => Promise<boolean>
  saveImageAs: (imagePath: string) => Promise<SaveImageResult>

  // Workspace
  getWorkspacePath: () => Promise<string>
}

declare global {
  interface Window {
    electronAPI: ElectronAPI
  }
}
