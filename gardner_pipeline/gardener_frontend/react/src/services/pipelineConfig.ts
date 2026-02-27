import { getLatestUserPreference } from './endpoints'
import type { RunPipelineRequest } from './types'

export interface PipelineSettings {
  skip_qc_filter: boolean
  skip_hvg: boolean
  skip_pca: boolean
  skip_neighbors: boolean
  min_genes: number
  min_cells: number
  pct_mt_max: number | null
  max_counts: number | null
  pct_hb_max: number | null
  n_top_genes_hvg: number
  flavor: 'seurat' | 'cell_ranger' | 'seurat_v3'
  target_sum: number
  n_comps: number
  svd_solver: 'arpack' | 'randomized'
  n_neighbors: number
  n_pcs: number
  clustering_method: 'leiden' | 'louvain' | 'cplearn'
  resolution: number
  run_hierarchical: boolean
  deg_groupby: string
  deg_method: 'wilcoxon' | 't-test' | 'logreg'
  n_top_genes_deg: number
  use_raw: boolean
  selected_annotation_models: string[]
  annotation_majority_voting: boolean
  annotation_top_n_genes: number
}

export interface PreferenceDefaults {
  settings: Partial<PipelineSettings>
  selectedAnnotationModels: string[]
  preferenceName?: string
  rawSettings: Record<string, unknown>
}

function compactDefined<T extends Record<string, unknown>>(obj: T): Partial<T> {
  return Object.fromEntries(
    Object.entries(obj).filter(([, value]) => value !== undefined)
  ) as Partial<T>
}

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

function normalizeModelSelections(settings: Record<string, unknown>): string[] {
  const flat = asStringArray(settings.selected_annotation_models)
  if (flat.length > 0) return flat

  const annotationParams = (settings.annotation_params ?? {}) as Record<string, unknown>
  const modelNames = asStringArray(annotationParams.model_names)
  const categories = asStringArray(annotationParams.categories)
  return [...modelNames, ...categories]
}

function splitAnnotationSelections(selected: string[]) {
  const normalized = selected
    .map((s) => (typeof s === 'string' ? s.trim() : String(s ?? '').trim()))
    .filter(Boolean)
  const modelNames = normalized.filter((s) => s.toLowerCase().endsWith('.pkl'))
  const categories = normalized.filter((s) => !s.toLowerCase().endsWith('.pkl'))
  return { normalized, modelNames, categories }
}

function fromFlatSettings(settings: Record<string, unknown>): Partial<PipelineSettings> {
  return {
    skip_qc_filter: asBoolean(settings.skip_qc_filter),
    skip_hvg: asBoolean(settings.skip_hvg),
    skip_pca: asBoolean(settings.skip_pca),
    skip_neighbors: asBoolean(settings.skip_neighbors),
    min_genes: asNumber(settings.min_genes),
    min_cells: asNumber(settings.min_cells),
    pct_mt_max: asNumber(settings.pct_mt_max) ?? null,
    max_counts: asNumber(settings.max_counts) ?? null,
    pct_hb_max: asNumber(settings.pct_hb_max) ?? null,
    n_top_genes_hvg: asNumber(settings.n_top_genes_hvg),
    flavor: (asString(settings.flavor) as PipelineSettings['flavor'] | undefined),
    target_sum: asNumber(settings.target_sum),
    n_comps: asNumber(settings.n_comps),
    svd_solver: (asString(settings.svd_solver) as PipelineSettings['svd_solver'] | undefined),
    n_neighbors: asNumber(settings.n_neighbors),
    n_pcs: asNumber(settings.n_pcs),
    clustering_method: (asString(settings.clustering_method) as PipelineSettings['clustering_method'] | undefined),
    resolution: asNumber(settings.resolution),
    run_hierarchical: asBoolean(settings.run_hierarchical),
    deg_groupby: asString(settings.deg_groupby),
    deg_method: (asString(settings.deg_method) as PipelineSettings['deg_method'] | undefined),
    n_top_genes_deg: asNumber(settings.n_top_genes_deg),
    use_raw: asBoolean(settings.use_raw),
    annotation_majority_voting: asBoolean(settings.annotation_majority_voting),
    annotation_top_n_genes: asNumber(settings.annotation_top_n_genes),
  }
}

function fromNestedParams(settings: Record<string, unknown>): Partial<PipelineSettings> {
  const preprocessing = (settings.preprocessing_params ?? {}) as Record<string, unknown>
  const clustering = (settings.clustering_params ?? {}) as Record<string, unknown>
  const dge = (settings.deg_params ?? {}) as Record<string, unknown>
  const annotation = (settings.annotation_params ?? {}) as Record<string, unknown>

  return {
    skip_qc_filter: asBoolean(preprocessing.skip_qc_filter),
    skip_hvg: asBoolean(preprocessing.skip_hvg),
    skip_pca: asBoolean(preprocessing.skip_pca),
    skip_neighbors: asBoolean(preprocessing.skip_neighbors),
    min_genes: asNumber(preprocessing.min_genes),
    min_cells: asNumber(preprocessing.min_cells),
    pct_mt_max: asNumber(preprocessing.pct_mt_max) ?? null,
    max_counts: asNumber(preprocessing.max_counts) ?? null,
    pct_hb_max: asNumber(preprocessing.pct_hb_max) ?? null,
    n_top_genes_hvg: asNumber(preprocessing.n_top_genes),
    flavor: (asString(preprocessing.flavor) as PipelineSettings['flavor'] | undefined),
    target_sum: asNumber(preprocessing.target_sum),
    n_comps: asNumber(preprocessing.n_comps),
    svd_solver: (asString(preprocessing.svd_solver) as PipelineSettings['svd_solver'] | undefined),
    n_neighbors: asNumber(preprocessing.n_neighbors),
    n_pcs: asNumber(preprocessing.n_pcs),
    clustering_method: (asString(clustering.method) as PipelineSettings['clustering_method'] | undefined),
    resolution: asNumber(clustering.resolution),
    run_hierarchical: asBoolean(clustering.run_hierarchical),
    deg_groupby: asString(dge.groupby),
    deg_method: (asString(dge.method) as PipelineSettings['deg_method'] | undefined),
    n_top_genes_deg: asNumber(dge.n_top_genes),
    use_raw: asBoolean(dge.use_raw),
    annotation_majority_voting: asBoolean(annotation.majority_voting),
    annotation_top_n_genes: asNumber(annotation.top_n_genes),
  }
}

export async function fetchLatestPreferenceDefaults(): Promise<PreferenceDefaults> {
  const pref = await getLatestUserPreference()
  const settings = pref.settings ?? {}
  const flat = compactDefined(fromFlatSettings(settings))
  const nested = compactDefined(fromNestedParams(settings))
  const merged: Partial<PipelineSettings> = {
    ...nested,
    ...flat,
  }

  return {
    settings: merged,
    selectedAnnotationModels: normalizeModelSelections(settings),
    preferenceName: pref.preference_name,
    rawSettings: settings,
  }
}

export function buildUpdatedPreferenceSettings(
  baseSettings: Record<string, unknown>,
  settings: Partial<PipelineSettings>,
  selectedAnnotationModels?: string[]
): Record<string, unknown> {
  const next = { ...baseSettings }
  const resolvedSelections = selectedAnnotationModels ?? normalizeModelSelections(baseSettings)
  const { normalized, modelNames, categories } = splitAnnotationSelections(resolvedSelections)
  const clusteringMethod =
    settings.clustering_method === 'leiden' ||
    settings.clustering_method === 'louvain' ||
    settings.clustering_method === 'cplearn'
      ? settings.clustering_method
      : undefined
  const targetClusterCol = clusteringMethod ?? 'leiden'

  const updatedFlat: Record<string, unknown> = {
    skip_qc_filter: settings.skip_qc_filter,
    skip_hvg: settings.skip_hvg,
    skip_pca: settings.skip_pca,
    skip_neighbors: settings.skip_neighbors,
    min_genes: settings.min_genes,
    min_cells: settings.min_cells,
    pct_mt_max: settings.pct_mt_max ?? null,
    max_counts: settings.max_counts ?? null,
    pct_hb_max: settings.pct_hb_max ?? null,
    n_top_genes_hvg: settings.n_top_genes_hvg,
    flavor: settings.flavor,
    target_sum: settings.target_sum,
    n_comps: settings.n_comps,
    svd_solver: settings.svd_solver,
    n_neighbors: settings.n_neighbors,
    n_pcs: settings.n_pcs,
    clustering_method: settings.clustering_method,
    resolution: settings.resolution,
    run_hierarchical: settings.run_hierarchical,
    deg_groupby: settings.deg_groupby,
    deg_method: settings.deg_method,
    n_top_genes_deg: settings.n_top_genes_deg,
    use_raw: settings.use_raw,
    annotation_majority_voting: settings.annotation_majority_voting,
    annotation_top_n_genes: settings.annotation_top_n_genes,
    selected_annotation_models: normalized,
  }

  for (const [k, v] of Object.entries(updatedFlat)) {
    if (v !== undefined) {
      next[k] = v
    }
  }

  const preprocessing = {
    ...((next.preprocessing_params as Record<string, unknown> | undefined) ?? {}),
    ...(settings.skip_qc_filter !== undefined ? { skip_qc_filter: settings.skip_qc_filter } : {}),
    ...(settings.skip_hvg !== undefined ? { skip_hvg: settings.skip_hvg } : {}),
    ...(settings.skip_pca !== undefined ? { skip_pca: settings.skip_pca } : {}),
    ...(settings.skip_neighbors !== undefined ? { skip_neighbors: settings.skip_neighbors } : {}),
    ...(settings.min_genes !== undefined ? { min_genes: settings.min_genes } : {}),
    ...(settings.min_cells !== undefined ? { min_cells: settings.min_cells } : {}),
    ...(settings.pct_mt_max !== undefined ? { pct_mt_max: settings.pct_mt_max ?? null } : {}),
    ...(settings.max_counts !== undefined ? { max_counts: settings.max_counts ?? null } : {}),
    ...(settings.pct_hb_max !== undefined ? { pct_hb_max: settings.pct_hb_max ?? null } : {}),
    ...(settings.n_top_genes_hvg !== undefined ? { n_top_genes: settings.n_top_genes_hvg } : {}),
    ...(settings.flavor !== undefined ? { flavor: settings.flavor } : {}),
    ...(settings.target_sum !== undefined ? { target_sum: settings.target_sum } : {}),
    ...(settings.n_comps !== undefined ? { n_comps: settings.n_comps } : {}),
    ...(settings.svd_solver !== undefined ? { svd_solver: settings.svd_solver } : {}),
    ...(settings.n_neighbors !== undefined ? { n_neighbors: settings.n_neighbors } : {}),
    ...(settings.n_pcs !== undefined ? { n_pcs: settings.n_pcs } : {}),
  }
  next.preprocessing_params = preprocessing

  const clustering = {
    ...((next.clustering_params as Record<string, unknown> | undefined) ?? {}),
    ...(settings.clustering_method !== undefined ? { method: settings.clustering_method } : {}),
    ...(settings.resolution !== undefined ? { resolution: settings.resolution } : {}),
    ...(settings.run_hierarchical !== undefined ? { run_hierarchical: settings.run_hierarchical } : {}),
  }
  next.clustering_params = clustering

  const deg = {
    ...((next.deg_params as Record<string, unknown> | undefined) ?? {}),
    ...(settings.deg_groupby !== undefined ? { groupby: settings.deg_groupby } : {}),
    ...(settings.deg_method !== undefined ? { method: settings.deg_method } : {}),
    ...(settings.n_top_genes_deg !== undefined ? { n_top_genes: settings.n_top_genes_deg } : {}),
    ...(settings.use_raw !== undefined ? { use_raw: settings.use_raw } : {}),
  }
  next.deg_params = deg

  const annotation = {
    ...((next.annotation_params as Record<string, unknown> | undefined) ?? {}),
    ...(settings.annotation_majority_voting !== undefined ? { majority_voting: settings.annotation_majority_voting } : {}),
    ...(settings.annotation_top_n_genes !== undefined ? { top_n_genes: settings.annotation_top_n_genes } : {}),
    ...(modelNames.length > 0 ? { model_names: modelNames } : { model_names: undefined }),
    ...(categories.length > 0 ? { categories } : { categories: undefined }),
    target_cluster_col: targetClusterCol,
  }
  next.annotation_params = annotation

  return next
}

export function settingsToRequest(
  projectId: string,
  datasetId: string,
  settings: PipelineSettings
): RunPipelineRequest {
  const inferredClusterCol =
    settings.clustering_method === 'leiden' ||
    settings.clustering_method === 'louvain' ||
    settings.clustering_method === 'cplearn'
      ? settings.clustering_method
      : 'leiden'
  const effectiveDegGroupby = settings.deg_groupby ?? inferredClusterCol

  const normalizedSelections = settings.selected_annotation_models
    .map((s) => s.trim())
    .filter(Boolean)

  const celltypistModels = normalizedSelections
    .filter((s) => s.toLowerCase().endsWith('.pkl'))

  const gseapyCategories = normalizedSelections
    .filter((s) => !s.toLowerCase().endsWith('.pkl'))

  return {
    project_id: projectId,
    dataset_id: datasetId,
    preprocessing_params: {
      skip_qc_filter: settings.skip_qc_filter,
      skip_hvg: settings.skip_hvg,
      skip_pca: settings.skip_pca,
      skip_neighbors: settings.skip_neighbors,
      min_genes: settings.min_genes,
      min_cells: settings.min_cells,
      pct_mt_max: settings.pct_mt_max ?? undefined,
      max_counts: settings.max_counts ?? undefined,
      pct_hb_max: settings.pct_hb_max ?? undefined,
      n_top_genes: settings.n_top_genes_hvg,
      flavor: settings.flavor,
      target_sum: settings.target_sum,
      n_comps: settings.n_comps,
      svd_solver: settings.svd_solver,
      n_neighbors: settings.n_neighbors,
      n_pcs: settings.n_pcs,
    },
    clustering_params: {
      method: settings.clustering_method,
      resolution: settings.resolution,
      run_hierarchical: settings.run_hierarchical,
    },
    deg_params: {
      groupby: effectiveDegGroupby,
      method: settings.deg_method,
      n_top_genes: settings.n_top_genes_deg,
    },
    annotation_params: {
      categories: gseapyCategories.length > 0 ? gseapyCategories : undefined,
      top_n_genes: settings.annotation_top_n_genes,
      model_names: celltypistModels.length > 0 ? celltypistModels : undefined,
      majority_voting: settings.annotation_majority_voting,
      target_cluster_col: inferredClusterCol,
    },
  }
}
