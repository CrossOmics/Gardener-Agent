import { useCallback, useEffect, useMemo, useState } from 'react'
import { Button, Form, Typography, message } from 'antd'
import { PlayCircleOutlined } from '@ant-design/icons'
import PipelineSettingsTabs from '../PipelineSetup/PipelineSettingsTabs'
import {
  getLatestStageSnapshot,
  runAnnotationFull,
  runClusteringCreate,
  runDGE,
  runPipeline,
  runPreprocessingFull,
} from '../../services/endpoints'
import { settingsToRequest, type PipelineSettings } from '../../services/pipelineConfig'
import { extractAnnotationModels, normalizeStageInputToFormValues } from '../FolderSummary/stagePrefill'
import './styles.css'

const { Title, Text } = Typography

interface DatasetStageRunnerProps {
  projectId: string
  datasetId: string
  datasetName: string
  onRunSuccess?: () => void
}

type PipelineTabKey = 'qc' | 'hvg' | 'pca' | 'neighbors' | 'clustering' | 'dge' | 'annotation'

function inferClusterMethodFromSnapshotId(snapshotId: string | undefined): 'leiden' | 'louvain' | 'cplearn' {
  const text = (snapshotId || '').toLowerCase()
  if (text.includes('louvain')) return 'louvain'
  if (text.includes('cplearn')) return 'cplearn'
  if (text.includes('leiden')) return 'leiden'
  return 'leiden'
}

export default function DatasetStageRunner({
  projectId,
  datasetId,
  datasetName,
  onRunSuccess,
}: DatasetStageRunnerProps) {
  const [form] = Form.useForm<PipelineSettings>()
  const [messageApi, contextHolder] = message.useMessage()
  const [activeTab, setActiveTab] = useState<PipelineTabKey>('qc')
  const [isRunningStage, setIsRunningStage] = useState(false)
  const [isRunningPipeline, setIsRunningPipeline] = useState(false)
  const [isLoadingDefaults, setIsLoadingDefaults] = useState(false)
  const [selectedAnnotationModels, setSelectedAnnotationModels] = useState<string[]>([])
  const [defaultFormValues, setDefaultFormValues] = useState<Partial<PipelineSettings>>({})
  const [clusteringSourceSnapshotId, setClusteringSourceSnapshotId] = useState<string>('')
  const [dgeSourceSnapshotId, setDgeSourceSnapshotId] = useState<string>('')
  const [annotationSourceSnapshotId, setAnnotationSourceSnapshotId] = useState<string>('')

  useEffect(() => {
    let cancelled = false
    const loadDefaults = async () => {
      setIsLoadingDefaults(true)
      const nextValues: Record<string, unknown> = {}
      let nextSelectedModels: string[] = []
      let preSnapshotId = ''
      let clusterSnapshotId = ''
      let dgeSnapshotId = ''

      try {
        const preprocessing = await getLatestStageSnapshot(datasetId, 'Preprocessing')
        if (!cancelled) {
          Object.assign(nextValues, normalizeStageInputToFormValues('preprocessing', preprocessing.input))
          preSnapshotId = preprocessing.snapshot_id
        }
      } catch {
        // ignore when absent
      }

      try {
        const clustering = await getLatestStageSnapshot(datasetId, 'Clustering')
        if (!cancelled) {
          const values = normalizeStageInputToFormValues('clustering', clustering.input)
          nextValues.clustering_method = values.method
          nextValues.resolution = values.resolution
          nextValues.run_hierarchical = values.run_hierarchical
          clusterSnapshotId = clustering.snapshot_id
          setClusteringSourceSnapshotId(clustering.parent_snapshot_id ?? preSnapshotId)
        }
      } catch {
        if (!cancelled) setClusteringSourceSnapshotId(preSnapshotId)
      }

      try {
        const dge = await getLatestStageSnapshot(datasetId, 'DGE')
        if (!cancelled) {
          const values = normalizeStageInputToFormValues('dge', dge.input)
          nextValues.deg_groupby = values.groupby
          nextValues.deg_method = values.method
          nextValues.n_top_genes_deg = values.n_top_genes
          nextValues.use_raw = values.use_raw
          dgeSnapshotId = dge.snapshot_id
          setDgeSourceSnapshotId(dge.parent_snapshot_id ?? clusterSnapshotId)
        }
      } catch {
        if (!cancelled) setDgeSourceSnapshotId(clusterSnapshotId)
      }

      try {
        const annotation = await getLatestStageSnapshot(datasetId, 'Annotation')
        if (!cancelled) {
          const values = normalizeStageInputToFormValues('annotation', annotation.input)
          nextValues.annotation_majority_voting = values.annotation_majority_voting
          nextValues.annotation_top_n_genes = values.annotation_top_n_genes
          const input = annotation.input && typeof annotation.input === 'object'
            ? (annotation.input as Record<string, unknown>)
            : {}
          nextSelectedModels = extractAnnotationModels(input)
          setAnnotationSourceSnapshotId(annotation.parent_snapshot_id ?? dgeSnapshotId)
        }
      } catch {
        if (!cancelled) setAnnotationSourceSnapshotId(dgeSnapshotId)
      }

      if (!cancelled) {
        const normalizedDefaults: Partial<PipelineSettings> = {
          ...(nextValues as Partial<PipelineSettings>),
        }
        form.setFieldsValue(normalizedDefaults)
        setDefaultFormValues(normalizedDefaults)
        setSelectedAnnotationModels(nextSelectedModels)
        setIsLoadingDefaults(false)
      }
    }

    void loadDefaults()
    return () => {
      cancelled = true
    }
  }, [datasetId, form])

  const preprocessingTabSet = useMemo(() => new Set<PipelineTabKey>(['qc', 'hvg', 'pca', 'neighbors']), [])

  const handleRunThisStage = useCallback(async () => {
    try {
      const values = await form.validateFields()
      const effectiveValues = {
        ...defaultFormValues,
        ...values,
      } as PipelineSettings
      setIsRunningStage(true)

      if (preprocessingTabSet.has(activeTab)) {
        const response = await runPreprocessingFull({
          project_id: projectId,
          dataset_id: datasetId,
          min_genes: effectiveValues.min_genes,
          min_cells: effectiveValues.min_cells,
          pct_mt_max: effectiveValues.pct_mt_max ?? undefined,
          cell_max_counts: effectiveValues.max_counts ?? undefined,
          pct_hb_max: effectiveValues.pct_hb_max ?? undefined,
          skip_qc_filter: effectiveValues.skip_qc_filter,
          n_top_genes: effectiveValues.n_top_genes_hvg,
          flavor: effectiveValues.flavor,
          target_sum: effectiveValues.target_sum,
          skip_hvg: effectiveValues.skip_hvg,
          n_comps: effectiveValues.n_comps,
          svd_solver: effectiveValues.svd_solver,
          skip_pca: effectiveValues.skip_pca,
          n_neighbors: effectiveValues.n_neighbors,
          n_pcs: effectiveValues.n_pcs,
          skip_neighbors: effectiveValues.skip_neighbors,
        })
        messageApi.success(response.message || response.msg || 'Preprocessing completed')
      } else if (activeTab === 'clustering') {
        if (!clusteringSourceSnapshotId) throw new Error('No source snapshot found for clustering')
        const response = await runClusteringCreate({
          project_id: projectId,
          dataset_id: datasetId,
          snapshot_id: clusteringSourceSnapshotId,
          method: effectiveValues.clustering_method,
          resolution: effectiveValues.resolution,
          run_hierarchical: effectiveValues.run_hierarchical,
        })
        messageApi.success(response.msg || 'Clustering completed')
      } else if (activeTab === 'dge') {
        if (!dgeSourceSnapshotId) throw new Error('No source snapshot found for DGE')
        const response = await runDGE({
          project_id: projectId,
          dataset_id: datasetId,
          snapshot_id: dgeSourceSnapshotId,
          groupby: effectiveValues.deg_groupby,
          method: effectiveValues.deg_method,
          n_top_genes: effectiveValues.n_top_genes_deg,
          use_raw: effectiveValues.use_raw,
        })
        messageApi.success(response.msg || 'DGE completed')
      } else if (activeTab === 'annotation') {
        if (!annotationSourceSnapshotId) throw new Error('No source snapshot found for annotation')
        const selected = selectedAnnotationModels.map((s) => s.trim()).filter(Boolean)
        const model_names = selected.filter((name) => name.toLowerCase().endsWith('.pkl'))
        const categories = selected.filter((name) => !name.toLowerCase().endsWith('.pkl'))
        const response = await runAnnotationFull({
          project_id: projectId,
          dataset_id: datasetId,
          snapshot_id: annotationSourceSnapshotId,
          majority_voting: effectiveValues.annotation_majority_voting,
          top_n_genes: effectiveValues.annotation_top_n_genes,
          model_names: model_names.length > 0 ? model_names : undefined,
          categories: categories.length > 0 ? categories : undefined,
          target_cluster_col: inferClusterMethodFromSnapshotId(dgeSourceSnapshotId),
        })
        messageApi.success(response.msg || 'Annotation completed')
      }

      onRunSuccess?.()
    } catch (error) {
      messageApi.error(error instanceof Error ? error.message : 'Run failed')
    } finally {
      setIsRunningStage(false)
    }
  }, [
    activeTab,
    annotationSourceSnapshotId,
    clusteringSourceSnapshotId,
    datasetId,
    dgeSourceSnapshotId,
    form,
    messageApi,
    onRunSuccess,
    preprocessingTabSet,
    projectId,
    selectedAnnotationModels,
    defaultFormValues,
  ])

  const handleRunWholePipeline = useCallback(async () => {
    try {
      const values = await form.validateFields()
      const effectiveValues = {
        ...defaultFormValues,
        ...values,
        selected_annotation_models: selectedAnnotationModels,
      } as PipelineSettings
      setIsRunningPipeline(true)
      const request = settingsToRequest(projectId, datasetId, effectiveValues)
      const response = await runPipeline(request)
      messageApi.success(response.msg || 'Pipeline completed successfully')
      onRunSuccess?.()
    } catch (error) {
      messageApi.error(error instanceof Error ? error.message : 'Pipeline execution failed')
    } finally {
      setIsRunningPipeline(false)
    }
  }, [datasetId, defaultFormValues, form, messageApi, onRunSuccess, projectId, selectedAnnotationModels])

  return (
    <div className="dataset-stage-runner">
      {contextHolder}
      <div className="dataset-stage-runner-header">
        <div className="dataset-stage-runner-title">
          <Title level={3}>Configure Analysis Pipeline</Title>
          <Text type="secondary">Dataset: {datasetName}</Text>
        </div>
        <div className="dataset-stage-runner-actions">
          <Button
            icon={<PlayCircleOutlined />}
            className="run-btn run-stage-btn"
            size="large"
            onClick={() => {
              void handleRunThisStage()
            }}
            loading={isRunningStage}
            disabled={isLoadingDefaults || isRunningPipeline}
          >
            Run This Stage
          </Button>
          <Button
            type="primary"
            icon={<PlayCircleOutlined />}
            className="run-btn"
            size="large"
            onClick={() => {
              void handleRunWholePipeline()
            }}
            loading={isRunningPipeline}
            disabled={isLoadingDefaults || isRunningStage}
          >
            Run Pipeline
          </Button>
        </div>
      </div>
      <div className="dataset-stage-runner-body">
        <Form form={form} layout="vertical">
          <PipelineSettingsTabs
            showAnnotationModelSelector
            selectedAnnotationModels={selectedAnnotationModels}
            onAnnotationSelectionChange={setSelectedAnnotationModels}
            onTabChange={(key) => setActiveTab(key as PipelineTabKey)}
          />
        </Form>
      </div>
    </div>
  )
}
