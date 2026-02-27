import { Form, InputNumber, Select, Switch, Tabs } from 'antd'
import type { FormInstance } from 'antd'
import { FormFieldLabel } from '../formShared'
import './styles.css'

interface PreprocessingSettingsTabsProps {
  form: FormInstance
  readOnly?: boolean
}

export default function PreprocessingSettingsTabs({ form, readOnly = false }: PreprocessingSettingsTabsProps) {
  const skipQcFilter = Boolean(Form.useWatch('skip_qc_filter', form))
  const skipHvg = Boolean(Form.useWatch('skip_hvg', form))
  const skipPca = Boolean(Form.useWatch('skip_pca', form))
  const skipNeighbors = Boolean(Form.useWatch('skip_neighbors', form))

  const tabItems = [
    {
      key: 'qc',
      label: 'QC Filtering',
      children: (
        <div className="settings-tab-content">
          <Form.Item
            name="skip_qc_filter"
            label={<FormFieldLabel text="Skip This Step" tooltip="Skip QC filtering and keep upstream data unchanged" />}
            valuePropName="checked"
          >
            <Switch disabled={readOnly} />
          </Form.Item>
          <Form.Item
            name="min_genes"
            label={<FormFieldLabel text="Minimum Genes per Cell" tooltip="Filter out cells with fewer genes than this threshold" />}
            rules={[{ required: !skipQcFilter }]}
          >
            <InputNumber className="settings-input" min={1} disabled={readOnly || skipQcFilter} />
          </Form.Item>
          <Form.Item
            name="min_cells"
            label={<FormFieldLabel text="Minimum Cells per Gene" tooltip="Filter out genes expressed in fewer cells than this threshold" />}
            rules={[{ required: !skipQcFilter }]}
          >
            <InputNumber className="settings-input" min={1} disabled={readOnly || skipQcFilter} />
          </Form.Item>
          <Form.Item
            name="pct_mt_max"
            label={<FormFieldLabel text="Maximum Mitochondrial %" tooltip="Maximum percentage of mitochondrial gene counts (filters dying cells)" />}
          >
            <InputNumber className="settings-input" min={0} max={100} placeholder="Optional" disabled={readOnly || skipQcFilter} />
          </Form.Item>
          <Form.Item
            name="max_counts"
            label={<FormFieldLabel text="Maximum Total Counts" tooltip="Maximum UMI counts per cell (filters potential doublets)" />}
          >
            <InputNumber className="settings-input" min={0} placeholder="Optional" disabled={readOnly || skipQcFilter} />
          </Form.Item>
          <Form.Item
            name="pct_hb_max"
            label={<FormFieldLabel text="Maximum Hemoglobin %" tooltip="Maximum percentage of hemoglobin gene counts" />}
          >
            <InputNumber className="settings-input" min={0} max={100} placeholder="Optional" disabled={readOnly || skipQcFilter} />
          </Form.Item>
        </div>
      ),
    },
    {
      key: 'hvg',
      label: 'HVG Selection',
      children: (
        <div className="settings-tab-content">
          <Form.Item
            name="skip_hvg"
            label={<FormFieldLabel text="Skip This Step" tooltip="Skip HVG selection and scaling" />}
            valuePropName="checked"
          >
            <Switch disabled={readOnly} />
          </Form.Item>
          <Form.Item
            name="n_top_genes_hvg"
            label={<FormFieldLabel text="Number of Top Genes" tooltip="Number of highly variable genes to select" />}
            rules={[{ required: !skipHvg }]}
          >
            <InputNumber className="settings-input" min={100} max={10000} disabled={readOnly || skipHvg} />
          </Form.Item>
          <Form.Item
            name="flavor"
            label={<FormFieldLabel text="Variability Method" tooltip="Algorithm for identifying highly variable genes" />}
            rules={[{ required: !skipHvg }]}
          >
            <Select
              disabled={readOnly || skipHvg}
              options={[
                { value: 'seurat', label: 'Seurat (Dispersion-based)' },
                { value: 'cell_ranger', label: 'Cell Ranger' },
                { value: 'seurat_v3', label: 'Seurat v3 (Requires raw counts)' },
              ]}
            />
          </Form.Item>
          <Form.Item
            name="target_sum"
            label={<FormFieldLabel text="Target Sum for Normalization" tooltip="Total counts per cell normalized to this value" />}
            rules={[{ required: !skipHvg }]}
          >
            <InputNumber className="settings-input" min={1} disabled={readOnly || skipHvg} />
          </Form.Item>
        </div>
      ),
    },
    {
      key: 'pca',
      label: 'PCA Analysis',
      children: (
        <div className="settings-tab-content">
          <Form.Item
            name="skip_pca"
            label={<FormFieldLabel text="Skip This Step" tooltip="Skip principal component analysis" />}
            valuePropName="checked"
          >
            <Switch disabled={readOnly} />
          </Form.Item>
          <Form.Item
            name="n_comps"
            label={<FormFieldLabel text="Number of Components" tooltip="Number of principal components to compute" />}
            rules={[{ required: !skipPca }]}
          >
            <InputNumber className="settings-input" min={1} max={200} disabled={readOnly || skipPca} />
          </Form.Item>
          <Form.Item
            name="svd_solver"
            label={<FormFieldLabel text="SVD Solver" tooltip="Algorithm for computing SVD" />}
            rules={[{ required: !skipPca }]}
          >
            <Select
              disabled={readOnly || skipPca}
              options={[
                { value: 'arpack', label: 'Arpack (Standard)' },
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
        <div className="settings-tab-content">
          <Form.Item
            name="skip_neighbors"
            label={<FormFieldLabel text="Skip This Step" tooltip="Skip neighbor graph construction" />}
            valuePropName="checked"
          >
            <Switch disabled={readOnly} />
          </Form.Item>
          <Form.Item
            name="n_neighbors"
            label={<FormFieldLabel text="Number of Neighbors" tooltip="Number of nearest neighbors for graph construction" />}
            rules={[{ required: !skipNeighbors }]}
          >
            <InputNumber className="settings-input" min={2} max={100} disabled={readOnly || skipNeighbors} />
          </Form.Item>
          <Form.Item
            name="n_pcs"
            label={<FormFieldLabel text="Number of Principal Components" tooltip="Number of PCs to use for neighbor graph" />}
            rules={[{ required: !skipNeighbors }]}
          >
            <InputNumber className="settings-input" min={1} max={100} disabled={readOnly || skipNeighbors} />
          </Form.Item>
        </div>
      ),
    },
  ]

  return <Tabs items={tabItems} size="small" />
}
