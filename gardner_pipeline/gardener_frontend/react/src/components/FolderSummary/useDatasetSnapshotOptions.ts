import { useEffect, useMemo, useState } from 'react'
import { getProjectDisplayMap } from '../../services/endpoints'
import type { DatasetEntityMap } from '../../services/types'
import type { SnapshotOption } from './formShared'

interface FileNode {
  name: string
  path: string
  isDirectory: boolean
  children?: FileNode[]
}

function normalizeFolderKey(name: string): string {
  return name.toLowerCase().replace(/[-\s]/g, '_')
}

function flattenDirectories(nodes: FileNode[]): FileNode[] {
  const result: FileNode[] = []
  for (const node of nodes) {
    if (!node.isDirectory) continue
    result.push(node)
    if (node.children && node.children.length > 0) {
      result.push(...flattenDirectories(node.children))
    }
  }
  return result
}

function getSnapshotOptions(datasetFiles: FileNode[], targetFolder: string): SnapshotOption[] {
  const dirs = flattenDirectories(datasetFiles)
  const pipelineDirs = dirs.filter((d) => normalizeFolderKey(d.name) === targetFolder)
  const snapshots = pipelineDirs.flatMap((d) => (d.children || []).filter((c) => c.isDirectory))
  return snapshots.map((s) => ({ label: s.name, value: s.name }))
}

function normalizeStageKey(name: string): string {
  return name.trim().toLowerCase()
}

function getSnapshotOptionsFromDisplayMap(
  datasetDisplay: DatasetEntityMap | null,
  stageName: 'Preprocessing' | 'Clustering' | 'DGE'
): SnapshotOption[] {
  if (!datasetDisplay) return []
  const target = normalizeStageKey(stageName)
  const stageEntry = Object.entries(datasetDisplay.snapshots).find(([key]) => normalizeStageKey(key) === target)
  if (!stageEntry) return []
  const [, snapshotMap] = stageEntry
  return Object.entries(snapshotMap)
    .filter(([, snapshotName]) => typeof snapshotName === 'string' && snapshotName.trim().length > 0)
    .map(([snapshotId, snapshotName]) => ({
      value: snapshotId,
      label: snapshotName,
    }))
}

export function useDatasetSnapshotOptions(projectId: string, datasetId: string) {
  const [datasetFiles, setDatasetFiles] = useState<FileNode[]>([])
  const [datasetDisplay, setDatasetDisplay] = useState<DatasetEntityMap | null>(null)

  useEffect(() => {
    const loadDatasetFiles = async () => {
      const electronApi = (
        window as unknown as {
          electronAPI?: {
            getDatasetFiles?: (projectId: string, datasetId: string) => Promise<FileNode[]>
          }
        }
      ).electronAPI

      if (!electronApi?.getDatasetFiles) {
        setDatasetFiles([])
        return
      }

      try {
        const [files, displayMap] = await Promise.all([
          electronApi.getDatasetFiles(projectId, datasetId),
          getProjectDisplayMap(projectId).catch(() => null),
        ])
        setDatasetFiles(files || [])
        setDatasetDisplay(displayMap?.datasets?.[datasetId] ?? null)
      } catch {
        setDatasetFiles([])
        setDatasetDisplay(null)
      }
    }

    loadDatasetFiles()
  }, [projectId, datasetId])

  const preprocessingSnapshotOptions = useMemo(
    () => {
      const fromMap = getSnapshotOptionsFromDisplayMap(datasetDisplay, 'Preprocessing')
      if (datasetDisplay) return fromMap
      return getSnapshotOptions(datasetFiles, 'preprocessing')
    },
    [datasetDisplay, datasetFiles]
  )
  const clusteringSnapshotOptions = useMemo(
    () => {
      const fromMap = getSnapshotOptionsFromDisplayMap(datasetDisplay, 'Clustering')
      if (datasetDisplay) return fromMap
      return getSnapshotOptions(datasetFiles, 'clustering')
    },
    [datasetDisplay, datasetFiles]
  )
  const dgeSnapshotOptions = useMemo(
    () => {
      const fromMap = getSnapshotOptionsFromDisplayMap(datasetDisplay, 'DGE')
      if (datasetDisplay) return fromMap
      return getSnapshotOptions(datasetFiles, 'dge')
    },
    [datasetDisplay, datasetFiles]
  )

  return {
    preprocessingSnapshotOptions,
    clusteringSnapshotOptions,
    dgeSnapshotOptions,
  }
}
