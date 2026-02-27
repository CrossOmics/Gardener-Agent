export type SettingsCategory = 'preprocessing' | 'clustering' | 'dge' | 'annotation'

const folderToCategoryMap: Record<string, SettingsCategory> = {
  preprocessing: 'preprocessing',
  preprocess: 'preprocessing',
  clustering: 'clustering',
  cluster_analysis: 'clustering',
  dge: 'dge',
  dge_analysis: 'dge',
  differential_expression: 'dge',
  annotation: 'annotation',
  cell_annotation: 'annotation',
}

const preprocessingSubfolders = new Set([
  'qc',
  'qc_filtering',
  'quality_control',
  'hvg',
  'hvg_selection',
  'highly_variable_genes',
  'pca',
  'pca_analysis',
  'neighbors',
  'neighbors_graph',
])

export function getCategoryFromFolderName(folderName: string): SettingsCategory | null {
  const normalized = folderName.toLowerCase().replace(/[-\s]/g, '_')
  return folderToCategoryMap[normalized] || null
}

export function isPreprocessingSubfolder(folderName: string): boolean {
  const normalized = folderName.toLowerCase().replace(/[-\s]/g, '_')
  return preprocessingSubfolders.has(normalized)
}
