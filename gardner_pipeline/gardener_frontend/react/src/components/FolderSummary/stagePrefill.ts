import type { SettingsCategory } from './folderCategory'
import type { SnapshotBranchName } from '../../services/types'

type AnyRecord = Record<string, unknown>

function asNumber(value: unknown): number | undefined {
  if (typeof value === 'number' && Number.isFinite(value)) return value
  if (typeof value === 'string' && value.trim() !== '' && Number.isFinite(Number(value))) {
    return Number(value)
  }
  return undefined
}

function asBoolean(value: unknown): boolean | undefined {
  if (typeof value === 'boolean') return value
  if (value === 'true') return true
  if (value === 'false') return false
  return undefined
}

function asString(value: unknown): string | undefined {
  if (typeof value !== 'string') return undefined
  const trimmed = value.trim()
  return trimmed ? trimmed : undefined
}

function asStringArray(value: unknown): string[] {
  if (!Array.isArray(value)) return []
  return value
    .map((item) => (typeof item === 'string' ? item.trim() : String(item ?? '').trim()))
    .filter(Boolean)
}

export function toBranchName(category: SettingsCategory): SnapshotBranchName {
  switch (category) {
    case 'preprocessing':
      return 'Preprocessing'
    case 'clustering':
      return 'Clustering'
    case 'dge':
      return 'DGE'
    case 'annotation':
      return 'Annotation'
  }
}

export function extractAnnotationModels(input: AnyRecord): string[] {
  const modelNames = asStringArray(input.model_names)
  const categories = asStringArray(input.categories)
  return [...modelNames, ...categories]
}

export function normalizeStageInputToFormValues(
  category: SettingsCategory,
  inputRaw: unknown
): AnyRecord {
  const input = (inputRaw && typeof inputRaw === 'object' ? inputRaw : {}) as AnyRecord

  if (category === 'preprocessing') {
    return {
      organism: asString(input.organism),
      custom_prefixes: input.custom_prefixes ?? null,
      skip_qc_calculation: asBoolean(input.skip_qc_calculation),
      skip_qc_filter: asBoolean(input.skip_qc_filter),
      skip_hvg: asBoolean(input.skip_hvg),
      skip_pca: asBoolean(input.skip_pca),
      skip_neighbors: asBoolean(input.skip_neighbors),
      min_genes: asNumber(input.min_genes),
      min_cells: asNumber(input.min_cells),
      pct_mt_max: asNumber(input.pct_mt_max) ?? null,
      max_counts: asNumber(input.max_counts) ?? asNumber(input.cell_max_counts) ?? null,
      pct_hb_max: asNumber(input.pct_hb_max) ?? null,
      n_top_genes_hvg: asNumber(input.n_top_genes),
      n_top_genes: asNumber(input.n_top_genes),
      flavor: asString(input.flavor),
      target_sum: asNumber(input.target_sum),
      n_comps: asNumber(input.n_comps),
      svd_solver: asString(input.svd_solver),
      n_neighbors: asNumber(input.n_neighbors),
      n_pcs: asNumber(input.n_pcs),
    }
  }

  if (category === 'clustering') {
    return {
      snapshot_id: asString(input.snapshot_id) ?? '',
      method: asString(input.method),
      resolution: asNumber(input.resolution),
      run_hierarchical: asBoolean(input.run_hierarchical),
    }
  }

  if (category === 'dge') {
    return {
      snapshot_id: asString(input.snapshot_id) ?? '',
      groupby: asString(input.groupby) ?? '',
      method: asString(input.method),
      n_top_genes: asNumber(input.n_top_genes),
      use_raw: asBoolean(input.use_raw),
    }
  }

  return {
    snapshot_id: asString(input.snapshot_id) ?? '',
    annotation_majority_voting: asBoolean(input.majority_voting) ?? true,
    annotation_top_n_genes: asNumber(input.top_n_genes) ?? 100,
  }
}
