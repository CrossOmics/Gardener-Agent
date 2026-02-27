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

export interface ProjectMetadata {
  name?: string
}
