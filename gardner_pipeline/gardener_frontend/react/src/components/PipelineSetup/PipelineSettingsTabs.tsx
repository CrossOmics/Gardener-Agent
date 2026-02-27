import { useState } from 'react'
import { Form, Input, InputNumber, Select, Switch, Tabs, Tooltip } from 'antd'
import { QuestionCircleOutlined } from '@ant-design/icons'
import AnnotationTools from '../AnnotationTools'

interface PipelineSettingsTabsProps {
  showAnnotationModelSelector?: boolean
  selectedAnnotationModels?: string[]
  onAnnotationSelectionChange?: (selectedNames: string[]) => void
  onTabChange?: (tabKey: string) => void
}

const Label = ({ text, tooltip }: { text: string; tooltip: string }) => (
  <span className="label-with-tooltip">
    {text}
    <Tooltip title={tooltip}>
      <QuestionCircleOutlined className="tooltip-icon" />
    </Tooltip>
  </span>
)

export default function PipelineSettingsTabs({
  showAnnotationModelSelector = false,
  selectedAnnotationModels = [],
  onAnnotationSelectionChange,
  onTabChange,
}: PipelineSettingsTabsProps) {
  const form = Form.useFormInstance()
  const [activeKey, setActiveKey] = useState('qc')
  const [annotationReloadToken, setAnnotationReloadToken] = useState(0)
  const skipQcFilter = Boolean(Form.useWatch('skip_qc_filter', form))
  const skipHvg = Boolean(Form.useWatch('skip_hvg', form))
  const skipPca = Boolean(Form.useWatch('skip_pca', form))
  const skipNeighbors = Boolean(Form.useWatch('skip_neighbors', form))

  const tabItems = [
    {
      key: 'qc',
      label: 'QC Filtering',
      children: (
        <div className="settings-form">
          <Form.Item
            name="skip_qc_filter"
            label={<Label text="Skip This Step" tooltip="Skip QC filtering and keep upstream data unchanged" />}
            valuePropName="checked"
          >
            <Switch />
          </Form.Item>
          <Form.Item
            name="min_genes"
            label={<Label text="Minimum Genes per Cell" tooltip="Minimum number of genes per cell (filters low-quality cells)" />}
            rules={[{ required: !skipQcFilter }, { type: 'number', min: 1 }]}
          >
            <InputNumber className="full-width-input" disabled={skipQcFilter} />
          </Form.Item>
          <Form.Item
            name="min_cells"
            label={<Label text="Minimum Cells per Gene" tooltip="Minimum number of cells per gene (filters low-expression genes)" />}
            rules={[{ required: !skipQcFilter }, { type: 'number', min: 1 }]}
          >
            <InputNumber className="full-width-input" disabled={skipQcFilter} />
          </Form.Item>
          <Form.Item
            name="pct_mt_max"
            label={<Label text="Maximum Mitochondrial %" tooltip="Maximum mitochondrial percentage (filters dead cells). Leave empty to skip." />}
          >
            <InputNumber className="full-width-input" min={0} max={100} placeholder="Optional" disabled={skipQcFilter} />
          </Form.Item>
          <Form.Item
            name="max_counts"
            label={<Label text="Maximum Total Counts" tooltip="Maximum total counts per cell (filters doublets). Leave empty to skip." />}
          >
            <InputNumber className="full-width-input" min={0} placeholder="Optional" disabled={skipQcFilter} />
          </Form.Item>
          <Form.Item
            name="pct_hb_max"
            label={<Label text="Maximum Hemoglobin %" tooltip="Maximum hemoglobin percentage (filters red blood cells). Leave empty to skip." />}
          >
            <InputNumber className="full-width-input" min={0} max={100} placeholder="Optional" disabled={skipQcFilter} />
          </Form.Item>
        </div>
      ),
    },
    {
      key: 'hvg',
      label: 'HVG Selection',
      children: (
        <div className="settings-form">
          <Form.Item
            name="skip_hvg"
            label={<Label text="Skip This Step" tooltip="Skip HVG selection and scaling" />}
            valuePropName="checked"
          >
            <Switch />
          </Form.Item>
          <Form.Item
            name="n_top_genes_hvg"
            label={<Label text="Number of Top Genes" tooltip="Number of highly variable genes to select" />}
            rules={[{ required: !skipHvg }, { type: 'number', min: 100 }]}
          >
            <InputNumber className="full-width-input" disabled={skipHvg} />
          </Form.Item>
          <Form.Item
            name="flavor"
            label={<Label text="Variability Method" tooltip="Method for computing gene variability" />}
            rules={[{ required: !skipHvg }]}
          >
            <Select
              disabled={skipHvg}
              options={[
                { value: 'seurat', label: 'Seurat (Dispersion-based, default)' },
                { value: 'cell_ranger', label: 'Cell Ranger' },
                { value: 'seurat_v3', label: 'Seurat v3 (Requires raw counts)' },
              ]}
            />
          </Form.Item>
          <Form.Item
            name="target_sum"
            label={<Label text="Target Sum for Normalization" tooltip="Target sum for normalization (total counts per cell normalized to this value)" />}
            rules={[{ required: !skipHvg }, { type: 'number', min: 1 }]}
          >
            <InputNumber className="full-width-input" disabled={skipHvg} />
          </Form.Item>
        </div>
      ),
    },
    {
      key: 'pca',
      label: 'PCA Analysis',
      children: (
        <div className="settings-form">
          <Form.Item
            name="skip_pca"
            label={<Label text="Skip This Step" tooltip="Skip principal component analysis" />}
            valuePropName="checked"
          >
            <Switch />
          </Form.Item>
          <Form.Item
            name="n_comps"
            label={<Label text="Number of Components" tooltip="Number of principal components" />}
            rules={[{ required: !skipPca }, { type: 'number', min: 1 }]}
          >
            <InputNumber className="full-width-input" disabled={skipPca} />
          </Form.Item>
          <Form.Item
            name="svd_solver"
            label={<Label text="SVD Solver" tooltip="SVD solver method" />}
            rules={[{ required: !skipPca }]}
          >
            <Select
              disabled={skipPca}
              options={[
                { value: 'arpack', label: 'Arpack (Standard, default)' },
                { value: 'randomized', label: 'Randomized (Faster for large datasets)' },
              ]}
            />
          </Form.Item>
        </div>
      ),
    },
    {
      key: 'neighbors',
      label: 'Neighbors Graph',
      children: (
        <div className="settings-form">
          <Form.Item
            name="skip_neighbors"
            label={<Label text="Skip This Step" tooltip="Skip neighbor graph construction" />}
            valuePropName="checked"
          >
            <Switch />
          </Form.Item>
          <Form.Item
            name="n_neighbors"
            label={<Label text="Number of Neighbors" tooltip="Number of neighbors for k-nearest neighbor graph" />}
            rules={[{ required: !skipNeighbors }, { type: 'number', min: 2 }]}
          >
            <InputNumber className="full-width-input" disabled={skipNeighbors} />
          </Form.Item>
          <Form.Item
            name="n_pcs"
            label={<Label text="Number of Principal Components" tooltip="Number of principal components used for building the neighbor graph" />}
            rules={[{ required: !skipNeighbors }, { type: 'number', min: 1 }]}
          >
            <InputNumber className="full-width-input" disabled={skipNeighbors} />
          </Form.Item>
        </div>
      ),
    },
    {
      key: 'clustering',
      label: 'Clustering',
      children: (
        <div className="settings-form">
          <Form.Item
            name="clustering_method"
            label={<Label text="Clustering Method" tooltip="Clustering algorithm" />}
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
            label={<Label text="Resolution" tooltip="Clustering resolution. Lower (0.1-0.3) = larger clusters, Higher (0.7-1.0) = smaller clusters" />}
            rules={[{ required: true }, { type: 'number', min: 0 }]}
          >
            <InputNumber className="full-width-input" step={0.1} />
          </Form.Item>
          <Form.Item
            name="run_hierarchical"
            label={<Label text="Run Hierarchical Clustering" tooltip="Whether to run hierarchical clustering analysis" />}
            valuePropName="checked"
          >
            <Switch />
          </Form.Item>
        </div>
      ),
    },
    {
      key: 'dge',
      label: 'DGE Analysis',
      children: (
        <div className="settings-form">
          <Form.Item
            name="deg_groupby"
            label={<Label text="Group By Column" tooltip="Grouping column used for differential expression. Default follows Clustering Method but can be edited freely." />}
          >
            <Input className="full-width-input" placeholder="e.g. leiden / louvain / custom_cluster_col" />
          </Form.Item>
          <Form.Item
            name="deg_method"
            label={<Label text="Statistical Test Method" tooltip="Statistical test method for identifying differentially expressed genes" />}
            rules={[{ required: true }]}
          >
            <Select
              options={[
                { value: 'wilcoxon', label: 'Wilcoxon (Rank-sum test, robust)' },
                { value: 't-test', label: "Student's t-test" },
                { value: 'logreg', label: 'Logistic Regression' },
              ]}
            />
          </Form.Item>
          <Form.Item
            name="n_top_genes_deg"
            label={<Label text="Number of Top Genes" tooltip="Number of top differentially expressed genes to visualize" />}
            rules={[{ required: true }, { type: 'number', min: 1 }]}
          >
            <InputNumber className="full-width-input" />
          </Form.Item>
          <Form.Item
            name="use_raw"
            label={<Label text="Use Raw Data" tooltip="Whether to use adata.raw (raw count data) for DEG calculation" />}
            valuePropName="checked"
          >
            <Switch />
          </Form.Item>
        </div>
      ),
    },
    {
      key: 'annotation',
      label: 'Annotation',
      children: showAnnotationModelSelector ? (
        <div className="annotation-tab-content">
          <div className="annotation-inline-row">
            <Form.Item
              name="annotation_majority_voting"
              className="annotation-voting-item"
              label={<Label text="Majority Voting" tooltip="Enable majority voting correction for CellTypist" />}
              valuePropName="checked"
            >
              <Switch />
            </Form.Item>
            <Form.Item
              name="annotation_top_n_genes"
              className="annotation-topn-item"
              label={<Label text="Top N Genes" tooltip="Number of top genes used for annotation (GSEApy)" />}
              rules={[{ required: true }, { type: 'number', min: 5, max: 500 }]}
            >
              <InputNumber />
            </Form.Item>
          </div>

          <h3 className="annotation-table-title">Model Selection</h3>

          <AnnotationTools
            embedded
            initialSelection={selectedAnnotationModels}
            reloadToken={annotationReloadToken}
            pageSize={5}
            onSelectionChange={onAnnotationSelectionChange}
          />
        </div>
      ) : (
        <div className="settings-form">
          <Form.Item
            name="annotation_majority_voting"
            label={<Label text="Majority Voting" tooltip="Enable majority voting correction for CellTypist" />}
            valuePropName="checked"
          >
            <Switch />
          </Form.Item>
          <Form.Item
            name="annotation_top_n_genes"
            label={<Label text="Top N Genes" tooltip="Number of top genes used for annotation (GSEApy)" />}
            rules={[{ required: true }, { type: 'number', min: 5, max: 500 }]}
          >
            <InputNumber className="full-width-input" />
          </Form.Item>
        </div>
      ),
    },
  ]

  return (
    <Tabs
      activeKey={activeKey}
      onChange={(nextKey) => {
        setActiveKey(nextKey)
        onTabChange?.(nextKey)
        if (nextKey === 'annotation') {
          setAnnotationReloadToken((prev) => prev + 1)
        }
      }}
      items={tabItems}
    />
  )
}
