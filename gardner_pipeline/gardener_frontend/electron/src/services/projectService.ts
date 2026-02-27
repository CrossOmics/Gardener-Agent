import fs from 'fs'
import path from 'path'
import { API_BASE } from '../config'
import type { DatasetInfo, FileNode, ProjectInfo, ProjectMetadata } from '../types'

function getProjectMetaPath(projectPath: string): string {
  return path.join(projectPath, 'project.meta.json')
}

function readProjectName(projectPath: string, fallbackId: string): string {
  try {
    const metadataPath = getProjectMetaPath(projectPath)
    if (!fs.existsSync(metadataPath)) {
      return fallbackId
    }
    const raw = fs.readFileSync(metadataPath, 'utf-8')
    const parsed = JSON.parse(raw) as ProjectMetadata
    const name = parsed?.name?.trim()
    return name || fallbackId
  } catch {
    return fallbackId
  }
}

function writeProjectName(projectPath: string, projectName: string) {
  const metadataPath = getProjectMetaPath(projectPath)
  const payload: ProjectMetadata = { name: projectName }
  fs.writeFileSync(metadataPath, JSON.stringify(payload, null, 2), 'utf-8')
}

function getDirectoryContentsWithPngOnly(dirPath: string, hiddenFolders: string[]): FileNode[] {
  try {
    const items = fs.readdirSync(dirPath)
    const result: FileNode[] = []

    for (const item of items) {
      if (hiddenFolders.includes(item)) {
        continue
      }

      const itemPath = path.join(dirPath, item)
      const stats = fs.statSync(itemPath)

      if (stats.isDirectory()) {
        const children = getDirectoryContentsWithPngOnly(itemPath, hiddenFolders)
        if (children.length > 0) {
          result.push({
            id: item,
            name: item,
            path: itemPath,
            isDirectory: true,
            children,
          })
        }
      } else if (item.toLowerCase().endsWith('.png')) {
        result.push({
          id: item,
          name: item,
          path: itemPath,
          isDirectory: false,
        })
      }
    }

    return result
  } catch (error) {
    console.error('Error reading directory:', error)
    return []
  }
}

export function createProjectService(projectsPath: string, hiddenFolders: string[]) {
  const getProjects = (): ProjectInfo[] => {
    try {
      if (!fs.existsSync(projectsPath)) {
        fs.mkdirSync(projectsPath, { recursive: true })
        return []
      }
      const items = fs.readdirSync(projectsPath)
      return items
        .filter(item => {
          const itemPath = path.join(projectsPath, item)
          return item.startsWith('p_') && fs.statSync(itemPath).isDirectory()
        })
        .map(item => {
          const itemPath = path.join(projectsPath, item)
          const stats = fs.statSync(itemPath)
          return {
            id: item,
            name: readProjectName(itemPath, item),
            path: itemPath,
            createdAt: stats.birthtime,
          }
        })
        .sort((a, b) => b.createdAt.getTime() - a.createdAt.getTime())
    } catch (error) {
      console.error('Error getting projects:', error)
      return []
    }
  }

  const getDatasets = (projectId: string): DatasetInfo[] => {
    try {
      const projectPath = path.join(projectsPath, projectId)
      if (!fs.existsSync(projectPath)) {
        return []
      }
      const items = fs.readdirSync(projectPath)
      return items
        .filter(item => {
          const itemPath = path.join(projectPath, item)
          return item.startsWith('ds_') && fs.statSync(itemPath).isDirectory()
        })
        .map(item => {
          const itemPath = path.join(projectPath, item)
          const stats = fs.statSync(itemPath)
          return {
            id: item,
            name: item,
            path: itemPath,
            createdAt: stats.birthtime,
          }
        })
        .sort((a, b) => b.createdAt.getTime() - a.createdAt.getTime())
    } catch (error) {
      console.error('Error getting datasets:', error)
      return []
    }
  }

  const getDatasetFiles = (projectId: string, datasetId: string): FileNode[] => {
    try {
      const datasetPath = path.join(projectsPath, projectId, datasetId)
      if (!fs.existsSync(datasetPath)) {
        return []
      }
      return getDirectoryContentsWithPngOnly(datasetPath, hiddenFolders)
    } catch (error) {
      console.error('Error getting dataset files:', error)
      return []
    }
  }

  const createProjectViaApi = async (projectName?: string): Promise<ProjectInfo | null> => {
    try {
      const normalizedName = projectName?.trim()
      if (!normalizedName) {
        return null
      }

      const response = await fetch(`${API_BASE}/project/create`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ project_name: normalizedName }),
      })

      if (!response.ok) {
        let backendMessage = `HTTP ${response.status}`
        try {
          const error = await response.json() as { detail?: string; message?: string }
          backendMessage = error.detail || error.message || backendMessage
        } catch {
          // Keep fallback message
        }
        throw new Error(backendMessage)
      }

      const payload = await response.json() as { project_id?: string; project_name?: string }
      const projectId = (payload.project_id || '').trim()
      if (!projectId) {
        throw new Error('Backend did not return project_id')
      }

      const projectPath = path.join(projectsPath, projectId)
      if (!fs.existsSync(projectPath)) {
        fs.mkdirSync(projectPath, { recursive: true })
      }

      const persistedName = (payload.project_name || normalizedName).trim() || projectId
      writeProjectName(projectPath, persistedName)
      const stats = fs.statSync(projectPath)
      return {
        id: projectId,
        name: persistedName,
        path: projectPath,
        createdAt: stats.birthtime,
      }
    } catch (error) {
      console.error('Error creating project via API:', error)
      return null
    }
  }

  const getProjectInfo = (projectId: string): ProjectInfo | null => {
    const projects = getProjects()
    return projects.find(p => p.id === projectId) || null
  }

  return {
    getProjects,
    getProjectInfo,
    getDatasets,
    getDatasetFiles,
    createProjectViaApi,
    getWorkspacePath: () => projectsPath,
  }
}
