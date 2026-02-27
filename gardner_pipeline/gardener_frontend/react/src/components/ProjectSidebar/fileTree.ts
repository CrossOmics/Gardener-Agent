export interface FileNode {
  id: string
  name: string
  displayName?: string
  path: string
  isDirectory: boolean
  children?: FileNode[]
}

export interface NavItem {
  type: 'project' | 'dataset' | 'folder' | 'snapshot' | 'image'
  id: string
  path: string
  name: string
  snapshotName?: string
  projectId?: string
  datasetId?: string
  snapshotId?: string
  isSnapshotFolder?: boolean
  snapshotCategory?: 'preprocessing' | 'clustering' | 'dge' | 'annotation'
}

const FOLDER_ORDER: Record<string, number> = {
  preprocessing: 1,
  clustering: 2,
  dge: 3,
  annotation: 4,
}

const SNAPSHOT_ID_PATTERNS = [/^s_[^\\/]+$/i, /^annot_[^\\/]+$/i, /^cluster_sub_[^\\/]+$/i]

export function inferSnapshotIdFromName(name: string): string | undefined {
  return SNAPSHOT_ID_PATTERNS.some((pattern) => pattern.test(name)) ? name : undefined
}

export function inferSnapshotIdFromDirectoryPath(dirPath: string): string | undefined {
  const normalized = dirPath.replace(/\\/g, '/')
  const segments = normalized.split('/').filter(Boolean)
  const last = segments[segments.length - 1]
  return last ? inferSnapshotIdFromName(last) : undefined
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

export function sortFiles(files: FileNode[]): FileNode[] {
  return [...files].sort((a, b) => {
    const aName = a.name.toLowerCase()
    const bName = b.name.toLowerCase()
    const aOrder = FOLDER_ORDER[aName] ?? 999
    const bOrder = FOLDER_ORDER[bName] ?? 999

    if (aOrder !== bOrder) {
      return aOrder - bOrder
    }

    return aName.localeCompare(bName)
  })
}

export function inferSnapshotIdFromPath(filePath: string): string | undefined {
  const normalized = filePath.replace(/\\/g, '/')
  const segments = normalized.split('/').filter(Boolean)
  for (let i = segments.length - 2; i >= 0; i--) {
    const seg = segments[i]
    if (SNAPSHOT_ID_PATTERNS.some((pattern) => pattern.test(seg))) {
      return seg
    }
  }
  return undefined
}

export function collectNavItems(files: FileNode[], projectId: string, datasetId: string): NavItem[] {
  const items: NavItem[] = []
  for (const file of sortFiles(files)) {
    if (file.isDirectory) {
      const category = inferSnapshotCategoryFromPath(file.path)
      const isSnapshotFolder = hasDirectImageChild(file)
      const snapshotId = isSnapshotFolder
        ? (inferSnapshotIdFromDirectoryPath(file.path) ?? inferSnapshotIdFromName(file.name) ?? file.id)
        : undefined
      const displayName = file.displayName ?? file.name
      items.push({
        type: isSnapshotFolder ? 'snapshot' : 'folder',
        id: file.id,
        path: file.path,
        name: displayName,
        snapshotName: isSnapshotFolder ? displayName : undefined,
        projectId,
        datasetId,
        snapshotId,
        isSnapshotFolder,
        snapshotCategory: isSnapshotFolder ? category : undefined,
      })
      if (file.children) {
        items.push(...collectNavItems(file.children, projectId, datasetId))
      }
      continue
    }

    items.push({
      type: 'image',
      id: file.id,
      path: file.path,
      name: file.name,
      projectId,
      datasetId,
      snapshotId: inferSnapshotIdFromPath(file.path),
    })
  }
  return items
}
