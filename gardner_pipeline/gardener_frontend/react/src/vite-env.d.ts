/// <reference types="vite/client" />

interface Project {
  id: string
  name: string
  path: string
  createdAt: Date
}

interface FileNode {
  id: string
  name: string
  path: string
  isDirectory: boolean
  children?: FileNode[]
}

interface Dataset {
  id: string
  name: string
  path: string
  createdAt: Date
}

interface DialogResult {
  canceled: boolean
  filePaths: string[]
}

interface SaveImageResult {
  canceled: boolean
  filePath?: string
}

interface ElectronAPI {
  getProjects: () => Promise<Project[]>
  getProjectInfo: (id: string) => Promise<Project | null>
  createProject: (projectName?: string) => Promise<Project | null>
  openProjectWindow: (projectId: string) => Promise<boolean>
  getDatasets: (projectId: string) => Promise<Dataset[]>
  getDatasetFiles: (projectId: string, datasetId: string) => Promise<FileNode[]>
  showUploadDialog: () => Promise<DialogResult>
  copyImageToClipboard: (imagePath: string) => Promise<boolean>
  saveImageAs: (imagePath: string) => Promise<SaveImageResult>
  getWorkspacePath: () => Promise<string>
}

declare global {
  interface Window {
    electronAPI: ElectronAPI
  }
}

export {}
