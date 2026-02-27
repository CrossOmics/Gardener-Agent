import { useState, useEffect, useCallback } from 'react'
import { Form, Input, InputNumber, Select, Switch, Button, message } from 'antd'
import { PlayCircleOutlined } from '@ant-design/icons'
import type { SettingsCategory } from '../folderCategory'
import {
  getLatestStageSnapshot,
  runAnnotationFull,
  runClusteringCreate,
  runDGE,
  runPreprocessingFull,
} from '../../../services/endpoints'
import AnnotationTools from '../../AnnotationTools'
import PreprocessingSettingsTabs from '../PreprocessingSettingsTabs'
import { FormFieldLabel, getMethodTypeColor, snapshotFilter, snapshotSort } from '../formShared'
import { useDatasetSnapshotOptions } from '../useDatasetSnapshotOptions'
import { extractAnnotationModels, normalizeStageInputToFormValues, toBranchName } from '../stagePrefill'
import './styles.css'

interface FolderSettingsPanelProps {
  projectId: string
  datasetId: string
  folderPath: string
  folderName: string
  category: SettingsCategory
  onRunSuccess?: () => void
  readOnly?: boolean
  readOnlyValues?: Record<string, unknown>
  readOnlySelectedMethods?: Array<{ id: string; name: string; type: 'celltypist' | 'gseapy' }>
}

function inferClusterMethodFromSnapshotId(snapshotId: string | undefined): 'leiden' | 'louvain' | 'cplearn' {
  const text = (snapshotId || '').toLowerCase()
  if (text.includes('louvain')) return 'louvain'
  if (text.includes('cplearn')) return 'cplearn'
  if (text.includes('leiden')) return 'leiden'
  return 'leiden'
}

export default function FolderSettingsPanel({
  projectId,
  datasetId,
  folderPath,
  folderName,
  category,
  onRunSuccess,
  readOnly = false,
  readOnlyValues,
  readOnlySelectedMethods = [],
}: FolderSettingsPanelProps) {
  const [form] = Form.useForm()
  const [messageApi, contextHolder] = message.useMessage()
  const [isRunning, setIsRunning] = useState(false)
  const [isPrefilling, setIsPrefilling] = useState(false)
  const [selectedAnnotationModels, setSelectedAnnotationModels] = useState<string[]>([])
  const { preprocessingSnapshotOptions, clusteringSnapshotOptions, dgeSnapshotOptions } =
    useDatasetSnapshotOptions(projectId, datasetId)

  const sourceSnapshotOptions =
    category === 'clustering'
      ? preprocessingSnapshotOptions
      : category === 'dge'
        ? clusteringSnapshotOptions
        : category === 'annotation'
          ? dgeSnapshotOptions
          : []

  useEffect(() => {
    if (readOnly) {
      form.setFieldsValue(readOnlyValues || {})
      return
    }
    let cancelled = false

    const loadFormValues = async () => {
      setIsPrefilling(true)
      let nextValues: Record<string, unknown> = {}
      let nextModels: string[] = []

      try {
        const latestStage = await getLatestStageSnapshot(datasetId, toBranchName(category))
        if (cancelled) return
        const normalized = normalizeStageInputToFormValues(category, latestStage.input)
        nextValues = { ...nextValues, ...normalized }
        if (category !== 'preprocessing' && latestStage.parent_snapshot_id) {
          nextValues.snapshot_id = latestStage.parent_snapshot_id
        }
        if (category === 'annotation') {
          const selectedFromInput = extractAnnotationModels(
            latestStage.input && typeof latestStage.input === 'object'
              ? (latestStage.input as Record<string, unknown>)
              : {}
          )
          if (selectedFromInput.length > 0) {
            nextModels = selectedFromInput
          }
        }
      } catch {
        // Keep empty values when no history exists yet.
      } finally {
        if (!cancelled) {
          form.setFieldsValue(nextValues)
          if (category === 'annotation') {
            setSelectedAnnotationModels(nextModels)
          }
          setIsPrefilling(false)
        }
      }
    }

    void loadFormValues()
    return () => {
      cancelled = true
    }
  }, [projectId, datasetId, folderPath, category, form, readOnly, readOnlyValues])

  useEffect(() => {
    if (category !== 'annotation') return
    const current = form.getFieldValue('snapshot_id')
    if (!current && dgeSnapshotOptions.length > 0) {
      form.setFieldValue('snapshot_id', dgeSnapshotOptions[0].value)
    }
  }, [category, dgeSnapshotOptions, form])

  useEffect(() => {
    if (readOnly || category === 'preprocessing') return
    const current = String(form.getFieldValue('snapshot_id') ?? '').trim()
    if (!current) return
    const exists = sourceSnapshotOptions.some((option) => option.value === current)
    if (!exists) {
      form.setFieldValue('snapshot_id', undefined)
    }
  }, [category, form, readOnly, sourceSnapshotOptions])

  const dgeSourceSnapshotId = Form.useWatch('snapshot_id', form)

  useEffect(() => {
    if (category !== 'dge') return
    const currentGroupby = String(form.getFieldValue('groupby') ?? '').trim()
    if (!currentGroupby) {
      form.setFieldValue('groupby', inferClusterMethodFromSnapshotId(String(dgeSourceSnapshotId ?? '')))
    }
  }, [category, dgeSourceSnapshotId, form])

  const handleAnnotationSelectionChange = useCallback((selectedNames: string[]) => {
    setSelectedAnnotationModels(selectedNames)
  }, [])

  const handleRun = async () => {
    try {
      let values: Record<string, unknown>
      if (category === 'annotation') {
        const annotationValues = await form.validateFields()
        values = {
          ...annotationValues,
          selected_models: selectedAnnotationModels,
        }
      } else {
        values = await form.validateFields()
      }

      setIsRunning(true)
      let runSucceeded = false

      try {
        if (category === 'preprocessing') {
          const response = await runPreprocessingFull({
            project_id: projectId,
            dataset_id: datasetId,
            min_genes: values.min_genes as number | undefined,
            min_cells: values.min_cells as number | undefined,
            pct_mt_max: (values.pct_mt_max as number | null | undefined) ?? undefined,
            cell_max_counts: (values.max_counts as number | null | undefined) ?? undefined,
            pct_hb_max: (values.pct_hb_max as number | null | undefined) ?? undefined,
            skip_qc_filter: values.skip_qc_filter as boolean | undefined,
            n_top_genes: (values.n_top_genes_hvg as number | undefined) ?? (values.n_top_genes as number | undefined),
            flavor: values.flavor as 'seurat' | 'cell_ranger' | 'seurat_v3' | undefined,
            target_sum: values.target_sum as number | undefined,
            skip_hvg: values.skip_hvg as boolean | undefined,
            n_comps: values.n_comps as number | undefined,
            svd_solver: values.svd_solver as string | undefined,
            skip_pca: values.skip_pca as boolean | undefined,
            n_neighbors: values.n_neighbors as number | undefined,
            n_pcs: values.n_pcs as number | undefined,
            skip_neighbors: values.skip_neighbors as boolean | undefined,
          })
          messageApi.success(response.message || response.msg || 'Preprocessing completed')
          runSucceeded = true
        } else if (category === 'clustering') {
          const response = await runClusteringCreate({
            project_id: projectId,
            dataset_id: datasetId,
            snapshot_id: values.snapshot_id as string | undefined,
            method: values.method as 'leiden' | 'louvain' | 'cplearn' | undefined,
            resolution: values.resolution as number | undefined,
            run_hierarchical: values.run_hierarchical as boolean | undefined,
          })
          messageApi.success(response.msg || 'Clustering completed')
          runSucceeded = true
        } else if (category === 'dge') {
          const response = await runDGE({
            project_id: projectId,
            dataset_id: datasetId,
            snapshot_id: values.snapshot_id as string | undefined,
            groupby: values.groupby as string | undefined,
            method: values.method as 'wilcoxon' | 't-test' | 'logreg' | undefined,
            n_top_genes: values.n_top_genes as number | undefined,
            use_raw: values.use_raw as boolean | undefined,
          })
          messageApi.success(response.msg || 'DGE completed')
          runSucceeded = true
        } else if (category === 'annotation') {
          const selected = selectedAnnotationModels
            .map((item) => item.trim())
            .filter(Boolean)

          const model_names = selected.filter((name) => name.toLowerCase().endsWith('.pkl'))
          const categories = selected
            .filter((name) => !name.toLowerCase().endsWith('.pkl'))

          const response = await runAnnotationFull({
            project_id: projectId,
            dataset_id: datasetId,
            snapshot_id: values.snapshot_id as string,
            majority_voting: values.annotation_majority_voting as boolean | undefined,
            top_n_genes: values.annotation_top_n_genes as number | undefined,
            model_names: model_names.length > 0 ? model_names : undefined,
            categories: categories.length > 0 ? categories : undefined,
            target_cluster_col: inferClusterMethodFromSnapshotId(values.snapshot_id as string | undefined),
          })
          messageApi.success(response.msg || 'Annotation completed')
          runSucceeded = true
        } else {
          messageApi.success(`Settings saved for "${folderName}"`)
        }
        if (runSucceeded) {
          onRunSuccess?.()
        }
      } catch (error) {
        messageApi.error(error instanceof Error ? error.message : 'Run failed')
      } finally {
        setIsRunning(false)
      }
    } catch {
      messageApi.error('Please check your input values')
    }
  }

  const renderClusteringForm = () => (
    <>
      <Form.Item
        name="snapshot_id"
        label={<FormFieldLabel text="Source Snapshot" tooltip="Select a preprocessing snapshot as the upstream input" />}
        rules={[{ required: true, message: 'Please select a preprocessing snapshot' }]}
      >
        {readOnly ? (
          <Input />
        ) : (
          <Select
            showSearch
            placeholder="Search and select preprocessing snapshot"
            options={preprocessingSnapshotOptions}
            filterOption={snapshotFilter}
            filterSort={snapshotSort}
          />
        )}
      </Form.Item>
      <Form.Item
        name="method"
        label={<FormFieldLabel text="Clustering Method" tooltip="Algorithm for clustering cells" />}
        rules={[{ required: true }]}
      >
        <Select
          options={[
            { value: 'leiden', label: 'Leiden (Modern, recommended)' },
            { value: 'louvain', label: 'Louvain (Classic)' },
            { value: 'cplearn', label: 'CPLearn (Core-periphery learning)' },
          ]}
        />
      </Form.Item>
      <Form.Item
        name="resolution"
        label={<FormFieldLabel text="Resolution" tooltip="Higher values produce more clusters" />}
        rules={[{ required: true }]}
      >
        <InputNumber className="settings-input" step={0.1} min={0.1} max={2.0} />
      </Form.Item>
      <Form.Item
        name="run_hierarchical"
        label={<FormFieldLabel text="Run Hierarchical" tooltip="Compute dendrogram for cluster hierarchy" />}
        valuePropName="checked"
      >
        <Switch />
      </Form.Item>
    </>
  )

  const renderDGEForm = () => (
    <>
      <Form.Item
        name="snapshot_id"
        label={<FormFieldLabel text="Source Snapshot" tooltip="Select a clustering snapshot as the upstream input" />}
        rules={[{ required: true, message: 'Please select a clustering snapshot' }]}
      >
        {readOnly ? (
          <Input />
        ) : (
          <Select
            showSearch
            placeholder="Search and select clustering snapshot"
            options={clusteringSnapshotOptions}
            filterOption={snapshotFilter}
            filterSort={snapshotSort}
          />
        )}
      </Form.Item>
      <Form.Item
        name="groupby"
        label={<FormFieldLabel text="Group By Column" tooltip="Grouping column for DGE. Default inferred from snapshot name, but editable as free text." />}
      >
        <Input placeholder="e.g. leiden / louvain / custom_cluster_col" />
      </Form.Item>
      <Form.Item
        name="method"
        label={<FormFieldLabel text="Statistical Test Method" tooltip="Method for differential expression testing" />}
        rules={[{ required: true }]}
      >
        <Select
          options={[
            { value: 'wilcoxon', label: 'Wilcoxon Rank-Sum (robust)' },
            { value: 't-test', label: "Student's t-test" },
            { value: 'logreg', label: 'Logistic Regression' },
          ]}
        />
      </Form.Item>
      <Form.Item
        name="n_top_genes"
        label={<FormFieldLabel text="Number of Top Genes" tooltip="Number of top marker genes per cluster" />}
        rules={[{ required: true }]}
      >
        <InputNumber className="settings-input" min={5} max={500} />
      </Form.Item>
      <Form.Item
        name="use_raw"
        label={<FormFieldLabel text="Use Raw Data" tooltip="Use raw counts for DE analysis (recommended)" />}
        valuePropName="checked"
      >
        <Switch />
      </Form.Item>
    </>
  )

  const renderAnnotationForm = () => (
    <div className="annotation-embedded-container">
      <Form.Item
        name="snapshot_id"
        label={<FormFieldLabel text="Source Snapshot" tooltip="Select a DGE snapshot as the upstream input" />}
        rules={[{ required: true, message: 'Please select a DGE snapshot' }]}
      >
        {readOnly ? (
          <Input />
        ) : (
          <Select
            showSearch
            placeholder="Search and select DGE snapshot"
            options={dgeSnapshotOptions}
            filterOption={snapshotFilter}
            filterSort={snapshotSort}
          />
        )}
      </Form.Item>

      <div className="annotation-inline-row">
        <Form.Item
          name="annotation_majority_voting"
          className="annotation-voting-item"
          label={<FormFieldLabel text="Majority Voting" tooltip="Enable majority voting correction for CellTypist" />}
          valuePropName="checked"
        >
          <Switch />
        </Form.Item>
        <Form.Item
          name="annotation_top_n_genes"
          className="annotation-topn-item"
          label={<FormFieldLabel text="Top N Genes" tooltip="Number of top genes used for annotation (GSEApy)" />}
          rules={[{ required: true }, { type: 'number', min: 5, max: 500 }]}
        >
          <InputNumber min={5} max={500} />
        </Form.Item>
      </div>

      <h3 className="annotation-table-title">Model Selection</h3>

      {readOnly ? (
        <div className="selected-section">
          <div className="selected-header">
            <h3 className="selected-title">Selected Methods</h3>
            <span className="color-legend">
              <span className="legend-item">
                <span className="legend-dot celltypist-dot"></span>CellTypist
              </span>
              <span className="legend-item">
                <span className="legend-dot gseapy-dot"></span>GSEApy
              </span>
            </span>
          </div>
          <div className="selected-models">
            {readOnlySelectedMethods.length === 0 ? (
              <span className="no-selection">No methods selected</span>
            ) : (
              readOnlySelectedMethods.map((model) => (
                <div
                  key={model.id}
                  className="selected-tag"
                  style={{ backgroundColor: getMethodTypeColor(model.type) }}
                >
                  <span className="tag-name">{model.name}</span>
                </div>
              ))
            )}
          </div>
        </div>
      ) : (
        <AnnotationTools
          embedded
          initialSelection={selectedAnnotationModels}
          onSelectionChange={handleAnnotationSelectionChange}
        />
      )}
    </div>
  )

  const renderSettingsForm = () => {
    switch (category) {
      case 'preprocessing':
        return <PreprocessingSettingsTabs form={form} readOnly={readOnly} />
      case 'clustering':
        return renderClusteringForm()
      case 'dge':
        return renderDGEForm()
      case 'annotation':
        return renderAnnotationForm()
      default:
        return null
    }
  }

  const isRunnable = !readOnly && (category === 'preprocessing' || category === 'clustering' || category === 'dge' || category === 'annotation')
  const buttonText = isRunnable ? 'Run' : 'Save Settings'
  const buttonIcon = isRunnable ? <PlayCircleOutlined /> : null

  return (
    <div className="folder-settings-panel">
      <div className="folder-settings-form">
        {contextHolder}
        {category === 'annotation' ? (
          <Form form={form} layout="vertical" size="small" disabled={readOnly}>
            {renderAnnotationForm()}
          </Form>
        ) : (
          <Form form={form} layout="vertical" size="small" disabled={readOnly}>
            {renderSettingsForm()}
          </Form>
        )}
        {!readOnly && (
          <div className="folder-settings-actions">
            <Button
              type="primary"
              icon={buttonIcon}
              onClick={handleRun}
              loading={isRunning || isPrefilling}
            >
              {buttonText}
            </Button>
          </div>
        )}
      </div>
    </div>
  )
}
