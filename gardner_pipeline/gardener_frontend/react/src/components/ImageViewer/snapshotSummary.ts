type AnyRecord = Record<string, unknown>

function pad2(n: number): string {
  return String(n).padStart(2, '0')
}

function formatCommonDateTime(raw: string): string {
  if (!raw) return 'N/A'
  const d = new Date(raw)
  if (Number.isNaN(d.getTime())) return raw
  const yyyy = d.getFullYear()
  const mm = pad2(d.getMonth() + 1)
  const dd = pad2(d.getDate())
  const hh = pad2(d.getHours())
  const mi = pad2(d.getMinutes())
  const ss = pad2(d.getSeconds())
  return `${yyyy}-${mm}-${dd} ${hh}:${mi}:${ss}`
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
  return undefined
}

function valueText(value: unknown): string {
  if (value === null || value === undefined) return 'null'
  if (typeof value === 'string') return value
  if (typeof value === 'number' || typeof value === 'boolean') return String(value)
  return JSON.stringify(value)
}

function getNested(obj: AnyRecord, path: string[]): unknown {
  let cur: unknown = obj
  for (const key of path) {
    if (!cur || typeof cur !== 'object') return undefined
    cur = (cur as AnyRecord)[key]
  }
  return cur
}

function listToEnglish(items: string[]): string {
  if (items.length === 0) return ''
  if (items.length === 1) return items[0]
  if (items.length === 2) return `${items[0]} and ${items[1]}`
  return `${items.slice(0, -1).join(', ')}, and ${items[items.length - 1]}`
}

function rangeSentence(
  minValue: unknown,
  maxValue: unknown,
  _minName: string,
  _maxName: string,
  unitSuffix = ''
): string {
  const min = asNumber(minValue)
  const max = asNumber(maxValue)
  if (min !== undefined && max === undefined) {
    return `at least ${min}${unitSuffix}`
  }
  if (min !== undefined && max !== undefined) {
    return `between ${min} and ${max}${unitSuffix}`
  }
  if (min === undefined && max !== undefined) {
    return `at most ${max}${unitSuffix}`
  }
  return 'not constrained'
}

function customPrefixesRender(input: AnyRecord): string {
  const customPrefixes = input.custom_prefixes
  if (customPrefixes === null || customPrefixes === undefined) return 'not set'
  return valueText(customPrefixes)
}

function enabledStepsSentence(input: AnyRecord): string {
  const steps: string[] = []
  if (asBoolean(input.skip_qc_calculation) === false) steps.push('compute QC metrics')
  if (asBoolean(input.skip_qc_filter) === false) steps.push('apply QC-based filtering')
  if (asBoolean(input.skip_hvg) === false) steps.push('select HVGs')
  if (asBoolean(input.skip_pca) === false) steps.push('run PCA')
  if (asBoolean(input.skip_neighbors) === false) steps.push('build the neighbor graph')
  if (steps.length === 0) return 'did not run any preprocessing sub-steps due to configuration'
  return `included ${listToEnglish(steps)}`
}

function hbFilterClause(input: AnyRecord): string {
  const pctHb = asNumber(input.pct_hb_max)
  if (pctHb === undefined) {
    return ', and no hemoglobin fraction threshold was applied'
  }
  return `, and hemoglobin fraction was filtered above ${pctHb}%`
}

const OUTPUT_FALLBACK_PATHS: Record<string, string[][]> = {
  n_cells_raw: [['summary', 'n_cells_raw']],
  n_cells_final: [['summary', 'n_cells_final']],
  retention_rate: [['summary', 'retention_rate']],
  n_genes_final: [['summary', 'n_genes_final']],
  median_genes_per_cell: [['summary', 'median_genes_per_cell']],
  median_counts_per_cell: [['summary', 'median_counts_per_cell']],
  median_pct_mt: [['qc_details', 'median_pct_mt']],
  low_quality_cells_removed: [['qc_details', 'filter_stats', 'low_quality_cells_removed']],
  n_hvg_selected: [['analysis_params', 'n_hvg_selected']],
  pcs_used: [['analysis_params', 'pcs_used']],
}

function readOutput(output: AnyRecord, key: string): string {
  const paths = OUTPUT_FALLBACK_PATHS[key] ?? [['summary', key]]
  for (const path of paths) {
    const value = getNested(output, path)
    if (value !== undefined) return valueText(value)
  }
  return valueText(output[key])
}

export function buildPreprocessingSummary(createTime: string, inputRaw: unknown, outputRaw: unknown): string {
  const input = (inputRaw && typeof inputRaw === 'object' ? inputRaw : {}) as AnyRecord
  const output = (outputRaw && typeof outputRaw === 'object' ? outputRaw : {}) as AnyRecord

  const organism = valueText(input.organism)
  const minGenes = input.min_genes
  const maxGenes = input.max_genes
  const cellMinCounts = input.cell_min_counts
  const cellMaxCounts = input.cell_max_counts
  const minCells = input.min_cells
  const maxCells = input.max_cells
  const geneMinCounts = input.gene_min_counts
  const geneMaxCounts = input.gene_max_counts
  const pctMtMax = valueText(input.pct_mt_max)
  const nTopGenes = valueText(input.n_top_genes)
  const flavor = valueText(input.flavor)
  const targetSum = valueText(input.target_sum)
  const nComps = valueText(input.n_comps)
  const svdSolver = valueText(input.svd_solver)
  const nNeighbors = valueText(input.n_neighbors)
  const nPcs = valueText(input.n_pcs)

  const cellGeneRange = rangeSentence(minGenes, maxGenes, 'min_genes', 'max_genes')
  const cellCountRange = rangeSentence(cellMinCounts, cellMaxCounts, 'cell_min_counts', 'cell_max_counts')
  const geneCellRange = rangeSentence(minCells, maxCells, 'min_cells', 'max_cells', ' cells')
  const geneCountRange = rangeSentence(geneMinCounts, geneMaxCounts, 'gene_min_counts', 'gene_max_counts')

  const displayTime = formatCommonDateTime(createTime)

  return [
    `This Preprocessing snapshot ran on ${organism} data at ${displayTime}, with custom feature prefixes ${customPrefixesRender(input)}.`,
    `It ${enabledStepsSentence(input)}.`,
    `For cell-level QC, cells were constrained by detected genes ${cellGeneRange} and total counts ${cellCountRange}.`,
    `For gene-level filtering, genes were constrained by detection in cells ${geneCellRange} and total counts ${geneCountRange}.`,
    `Cells were filtered if mitochondrial fraction exceeded ${pctMtMax}%${hbFilterClause(input)}.`,
    `HVGs (highly variable genes, genes with the strongest variation across cells) were selected as ${nTopGenes} using ${flavor}, counts were normalized to ${targetSum} per cell, PCA computed ${nComps} PCs with ${svdSolver}, and the neighbor graph used ${nNeighbors} neighbors over ${nPcs} PCs.`,
    `Results retained ${readOutput(output, 'n_cells_final')}/${readOutput(output, 'n_cells_raw')} cells (${readOutput(output, 'retention_rate')}), with ${readOutput(output, 'n_genes_final')} genes, median ${readOutput(output, 'median_genes_per_cell')} genes/cell and ${readOutput(output, 'median_counts_per_cell')} counts/cell.`,
    `QC diagnostics (per-cell quality metrics such as detected genes, total counts, and mitochondrial fraction) reported median mitochondrial fraction ${readOutput(output, 'median_pct_mt')}% and removed ${readOutput(output, 'low_quality_cells_removed')} low-quality cells; downstream settings recorded ${readOutput(output, 'n_hvg_selected')} HVGs and ${readOutput(output, 'pcs_used')} PCs.`,
  ].join(' ')
}

export function buildClusteringSummary(createTime: string, inputRaw: unknown, outputRaw: unknown): string {
  const input = (inputRaw && typeof inputRaw === 'object' ? inputRaw : {}) as AnyRecord
  const output = (outputRaw && typeof outputRaw === 'object' ? outputRaw : {}) as AnyRecord
  const displayTime = formatCommonDateTime(createTime)

  const method = valueText(getNested(output, ['summary', 'method']) ?? input.method)
  const resolution = valueText(getNested(output, ['summary', 'resolution']) ?? input.resolution)
  const runHierarchical = asBoolean(input.run_hierarchical)
  const hierarchicalClause = runHierarchical
    ? ', with hierarchical clustering enabled'
    : ', with hierarchical clustering disabled'

  const nClusters = readOutput(output, 'n_clusters')
  const nCellsTotal = readOutput(output, 'n_cells_total')

  return `This Clustering snapshot ran at ${displayTime} using the ${method} algorithm with resolution ${resolution}${hierarchicalClause}. Results identified ${nClusters} clusters across ${nCellsTotal} cells.`
}
