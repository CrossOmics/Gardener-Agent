import { Tooltip } from 'antd'
import { QuestionCircleOutlined } from '@ant-design/icons'

export interface SnapshotOption {
  label: string
  value: string
}

export function FormFieldLabel({ text, tooltip }: { text: string; tooltip: string }) {
  return (
    <span style={{ display: 'inline-flex', alignItems: 'center', gap: 4 }}>
      {text}
      <Tooltip title={tooltip}>
        <QuestionCircleOutlined style={{ color: '#999', fontSize: 12 }} />
      </Tooltip>
    </span>
  )
}

export function snapshotFilter(input: string, option?: { label?: string; value?: string }) {
  const keyword = input.trim().toLowerCase()
  if (!keyword) return true
  const label = String(option?.label ?? '').toLowerCase()
  const value = String(option?.value ?? '').toLowerCase()
  return label.includes(keyword) || value.includes(keyword)
}

export function snapshotSort(
  a: { label?: string; value?: string },
  b: { label?: string; value?: string },
  info: { searchValue: string }
) {
  const keyword = info.searchValue.trim().toLowerCase()
  const aLabel = String(a.label ?? '').toLowerCase()
  const bLabel = String(b.label ?? '').toLowerCase()
  if (!keyword) return aLabel.localeCompare(bLabel)

  const aStarts = aLabel.startsWith(keyword) ? 1 : 0
  const bStarts = bLabel.startsWith(keyword) ? 1 : 0
  if (aStarts !== bStarts) return bStarts - aStarts

  const aIncludes = aLabel.includes(keyword) ? 1 : 0
  const bIncludes = bLabel.includes(keyword) ? 1 : 0
  if (aIncludes !== bIncludes) return bIncludes - aIncludes

  return aLabel.localeCompare(bLabel)
}

export function getMethodTypeColor(type: 'celltypist' | 'gseapy') {
  return type === 'celltypist' ? '#fa8c16' : '#1890ff'
}
