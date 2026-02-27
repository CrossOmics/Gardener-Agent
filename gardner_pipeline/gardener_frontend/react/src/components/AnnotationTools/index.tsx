import { useState, useMemo, useEffect } from 'react'
import { Input, Select, Table, Button, Checkbox, message } from 'antd'
import { SearchOutlined, CloseOutlined } from '@ant-design/icons'
import type { ColumnsType } from 'antd/es/table'
import { searchAnnotationOptions, updateUserPreference } from '../../services/endpoints'
import type { AnnotationSearchItem } from '../../services/types'
import { buildUpdatedPreferenceSettings, fetchLatestPreferenceDefaults } from '../../services/pipelineConfig'
import './styles.css'

type SearchFilter = 'all' | 'celltypist' | 'gseapy'
type ModelType = 'celltypist' | 'gseapy'

interface AnnotationModel {
  id: string
  name: string
  description: string
  type: ModelType
}

const DEFAULT_PAGE_SIZE = 7
const EMPTY_SELECTION: string[] = []

function normalizeSelection(values: string[] | undefined): string[] {
  if (!values) return []
  return values
    .map((item) => (typeof item === 'string' ? item.trim() : String(item ?? '').trim()))
    .filter(Boolean)
}

interface AnnotationToolsProps {
  /** When true, hides header and action buttons for embedding in other components */
  embedded?: boolean
  /** Optional project context for loading/saving defaults */
  projectId?: string
  /** Initial selected model names, mainly for embedded mode */
  initialSelection?: string[] | undefined
  /** Callback when selection changes */
  onSelectionChange?: (selectedNames: string[]) => void
  /** Changing token forces list re-fetch */
  reloadToken?: number
  /** Table page size */
  pageSize?: number
}

export default function AnnotationTools({
  embedded = false,
  projectId,
  initialSelection,
  onSelectionChange,
  reloadToken = 0,
  pageSize = DEFAULT_PAGE_SIZE,
}: AnnotationToolsProps = {}) {
  const [searchText, setSearchText] = useState('')
  const [searchFilter, setSearchFilter] = useState<SearchFilter>('all')
  const [selectedNames, setSelectedNames] = useState<Set<string>>(() => {
    if (initialSelection !== undefined) {
      return new Set(normalizeSelection(initialSelection))
    }
    return new Set(EMPTY_SELECTION)
  })
  const [models, setModels] = useState<AnnotationModel[]>([])
  const [isLoading, setIsLoading] = useState(false)
  const [isSaving, setIsSaving] = useState(false)
  const [preferenceName, setPreferenceName] = useState<string>('')
  const [baseSettings, setBaseSettings] = useState<Record<string, unknown>>({})
  const [modelRegistry, setModelRegistry] = useState<Record<string, AnnotationModel>>({})
  const [messageApi, contextHolder] = message.useMessage()

  // Sync with parent-controlled initial selection (avoid no-op loops)
  useEffect(() => {
    if (initialSelection === undefined) {
      return
    }
    const normalized = normalizeSelection(initialSelection)
    setSelectedNames((prev) => {
      const next = new Set(normalized)
      if (prev.size === next.size) {
        let same = true
        for (const id of prev) {
          if (!next.has(id)) {
            same = false
            break
          }
        }
        if (same) return prev
      }
      return next
    })
  }, [initialSelection])

  // Load project defaults in standalone mode
  useEffect(() => {
    if (embedded || !projectId || initialSelection !== undefined) {
      return
    }
    let mounted = true
    void (async () => {
      try {
        const defaults = await fetchLatestPreferenceDefaults()
        if (!mounted) return
        setSelectedNames(new Set(normalizeSelection(defaults.selectedAnnotationModels)))
        setPreferenceName(defaults.preferenceName ?? '')
        setBaseSettings(defaults.rawSettings)
      } catch (error) {
        if (!mounted) return
        console.error('Failed to load latest preferences for annotation defaults:', error)
      }
    })()
    return () => {
      mounted = false
    }
  }, [embedded, projectId, initialSelection])

  // Notify parent of selection changes
  useEffect(() => {
    onSelectionChange?.(Array.from(selectedNames))
  }, [selectedNames, onSelectionChange])

  const mapSearchItem = (item: AnnotationSearchItem): AnnotationModel => ({
    id: String(item.id ?? ''),
    name: String(item.method_name ?? '').trim(),
    description: typeof item.description === 'string' ? item.description : '',
    type: item.type === 'celltypist' ? 'celltypist' : 'gseapy',
  })

  const fetchModels = async (keyword: string, filter: SearchFilter) => {
    setIsLoading(true)
    try {
      const response = await searchAnnotationOptions({
        keyword: keyword.trim() || undefined,
        type: filter,
      })
      const mapped = response
        .map(mapSearchItem)
        .filter((item) => item.id !== '' && item.name !== '')
      setModels(mapped)
      setModelRegistry((prev) => {
        const next = { ...prev }
        for (const model of mapped) {
          next[model.id] = model
        }
        return next
      })
    } catch (error) {
      console.error('Failed to search annotation methods:', error)
      setModels([])
    } finally {
      setIsLoading(false)
    }
  }

  useEffect(() => {
    void fetchModels('', 'all')
  }, [reloadToken])

  const handleSearch = () => {
    void fetchModels(searchText, searchFilter)
  }

  const selectedLowerSet = useMemo(
    () => new Set(Array.from(selectedNames).map((name) => String(name).toLowerCase())),
    [selectedNames]
  )

  const selectedModels = useMemo(() => {
    const byName = Object.values(modelRegistry).reduce<Record<string, AnnotationModel>>((acc, model) => {
      acc[String(model.name).toLowerCase()] = model
      return acc
    }, {})
    return Array.from(selectedNames).map((name) => {
      const normalizedName = String(name)
      const matched = byName[normalizedName.toLowerCase()]
      if (matched) return matched
      return {
        id: `manual-${normalizedName}`,
        name: normalizedName,
        description: '',
        type: normalizedName.toLowerCase().endsWith('.pkl') ? 'celltypist' : 'gseapy',
      } as AnnotationModel
    })
  }, [selectedNames, modelRegistry])

  const handleToggleSelect = (name: string) => {
    const normalizedName = String(name ?? '').trim()
    if (!normalizedName) return
    setSelectedNames((prev) => {
      const next = new Set(prev)
      const existing = Array.from(next).find((item) => String(item).toLowerCase() === normalizedName.toLowerCase())
      if (existing) {
        next.delete(existing)
      } else {
        next.add(normalizedName)
      }
      return next
    })
  }

  const handleRemoveSelected = (name: string) => {
    const normalizedName = String(name ?? '').trim()
    if (!normalizedName) return
    setSelectedNames((prev) => {
      const next = new Set(prev)
      const existing = Array.from(next).find((item) => String(item).toLowerCase() === normalizedName.toLowerCase())
      if (existing) {
        next.delete(existing)
      }
      return next
    })
  }

  const handleSave = async () => {
    try {
      if (!preferenceName) {
        messageApi.error('Missing latest preference name, cannot update defaults')
        return
      }
      setIsSaving(true)
      const selected = Array.from(selectedNames)
      const settings = buildUpdatedPreferenceSettings(baseSettings, {}, selected)
      await updateUserPreference({
        name: preferenceName,
        settings,
      })
      setBaseSettings(settings)
      messageApi.success('Annotation defaults updated')
    } catch (error) {
      messageApi.error(error instanceof Error ? error.message : 'Failed to update annotation defaults')
    } finally {
      setIsSaving(false)
    }
  }

  const handleDiscard = () => {
    setSelectedNames(new Set())
    setSearchText('')
    setSearchFilter('all')
    void fetchModels('', 'all')
  }

  const getTypeColor = (type: ModelType) => {
    return type === 'celltypist' ? '#fa8c16' : '#1890ff'
  }

  const columns: ColumnsType<AnnotationModel> = [
    {
      title: 'Model Name',
      dataIndex: 'name',
      key: 'name',
      render: (name: string, record) => (
        <span className="model-name" style={{ color: getTypeColor(record.type) }}>
          {name}
        </span>
      ),
    },
    {
      title: 'Description',
      dataIndex: 'description',
      key: 'description',
      ellipsis: true,
    },
    {
      title: 'Selected',
      key: 'selected',
      width: 100,
      align: 'center',
      render: (_, record) => (
        <Checkbox
          checked={selectedLowerSet.has(record.name.toLowerCase())}
          onChange={() => handleToggleSelect(record.name)}
        />
      ),
    },
  ]

  return (
    <div className={`annotation-container ${embedded ? 'embedded' : ''}`}>
      {contextHolder}
      {!embedded && (
        <div className="annotation-header">
          <h2 className="annotation-title">Annotation Tools</h2>
        </div>
      )}

      <div className="annotation-content">
        <div className="search-bar">
          <Input
            placeholder="Search models..."
            addonBefore={
              <Select
                value={searchFilter}
                onChange={setSearchFilter}
                options={[
                  { value: 'all', label: 'All Types' },
                  { value: 'celltypist', label: 'CellTypist' },
                  { value: 'gseapy', label: 'GSEApy' },
                ]}
                className="search-filter-select"
              />
            }
            value={searchText}
            onChange={(e) => setSearchText(e.target.value)}
            onPressEnter={handleSearch}
            className="search-input"
            allowClear
          />
          <Button
            type="primary"
            shape="circle"
            icon={<SearchOutlined />}
            onClick={handleSearch}
            className="search-btn"
          />
        </div>

        <Table
          columns={columns}
          dataSource={models}
          rowKey="id"
          loading={isLoading}
          pagination={{
            pageSize,
            showSizeChanger: false,
            showTotal: (total) => `Total ${total} models`,
          }}
          className="models-table"
          rowClassName={(record) =>
            record.type === 'celltypist' ? 'row-celltypist' : 'row-gseapy'
          }
        />

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
            {selectedModels.length === 0 ? (
              <span className="no-selection">No methods selected</span>
            ) : (
              selectedModels.map((model) => (
                <div
                  key={model.id}
                  className="selected-tag"
                  style={{ backgroundColor: getTypeColor(model.type) }}
                >
                  <span className="tag-name">{model.name}</span>
                  <CloseOutlined
                    className="tag-close"
                    onClick={() => handleRemoveSelected(model.name)}
                  />
                </div>
              ))
            )}
          </div>
        </div>

        {!embedded && (
          <div className="action-buttons">
            <Button
              type="primary"
              className="save-btn"
              onClick={handleSave}
              disabled={selectedNames.size === 0 || isSaving}
              loading={isSaving}
            >
              Save
            </Button>
            <Button className="discard-btn" onClick={handleDiscard}>
              Discard
            </Button>
          </div>
        )}
      </div>
    </div>
  )
}
