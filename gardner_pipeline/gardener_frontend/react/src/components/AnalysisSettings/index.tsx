import { useEffect, useRef, useState } from 'react'
import { Form, Button, message } from 'antd'
import { SaveOutlined } from '@ant-design/icons'
import {
  buildUpdatedPreferenceSettings,
  type PipelineSettings,
  fetchLatestPreferenceDefaults,
} from '../../services/pipelineConfig'
import { updateUserPreference } from '../../services/endpoints'
import PipelineSettingsTabs from '../PipelineSetup/PipelineSettingsTabs'
import './styles.css'

type ProjectDefaultSettings = Omit<PipelineSettings, 'selected_annotation_models'>

function inferGroupByFromMethod(method: PipelineSettings['clustering_method'] | undefined): string {
  if (method === 'louvain' || method === 'cplearn' || method === 'leiden') return method
  return 'leiden'
}

interface AnalysisSettingsProps {
  projectId: string
}

export default function AnalysisSettings({ projectId }: AnalysisSettingsProps) {
  const [form] = Form.useForm<ProjectDefaultSettings>()
  const [messageApi, contextHolder] = message.useMessage()
  const [isSaving, setIsSaving] = useState(false)
  const [preferenceName, setPreferenceName] = useState<string>('')
  const [baseSettings, setBaseSettings] = useState<Record<string, unknown>>({})
  const clusteringMethod = Form.useWatch('clustering_method', form)
  const prevInferredGroupByRef = useRef<string>('leiden')

  useEffect(() => {
    let mounted = true
    void (async () => {
      try {
        const { settings, preferenceName: loadedPreferenceName, rawSettings } = await fetchLatestPreferenceDefaults()
        if (!mounted) return
        form.setFieldsValue(settings)
        setPreferenceName(loadedPreferenceName ?? '')
        setBaseSettings(rawSettings)
        prevInferredGroupByRef.current = inferGroupByFromMethod(settings.clustering_method)
      } catch (error) {
        if (!mounted) return
        const msg = error instanceof Error ? error.message : 'Failed to load latest preferences'
        messageApi.error(msg)
      }
    })()
    return () => {
      mounted = false
    }
  }, [projectId, form])

  useEffect(() => {
    const nextInferred = inferGroupByFromMethod(clusteringMethod)
    const current = form.getFieldValue('deg_groupby')
    if (current === undefined || current === '' || current === prevInferredGroupByRef.current) {
      form.setFieldValue('deg_groupby', nextInferred)
    }
    prevInferredGroupByRef.current = nextInferred
  }, [clusteringMethod, form])

  const handleSave = async () => {
    try {
      const values = await form.validateFields()
      if (!preferenceName) {
        messageApi.error('Missing latest preference name, cannot update defaults')
        return
      }
      setIsSaving(true)
      const settings = buildUpdatedPreferenceSettings(baseSettings, values as Partial<PipelineSettings>)
      await updateUserPreference({
        name: preferenceName,
        settings,
      })
      setBaseSettings(settings)
      messageApi.success('Settings saved successfully')
    } catch (error) {
      if (error instanceof Error) {
        messageApi.error(error.message)
      } else {
        messageApi.error('Please check your input values')
      }
    } finally {
      setIsSaving(false)
    }
  }

  return (
    <div className="settings-container">
      {contextHolder}
      <div className="settings-header">
        <h2 className="settings-title">Customize Analysis</h2>
        <Button type="primary" icon={<SaveOutlined />} onClick={handleSave} loading={isSaving}>
          Save Settings
        </Button>
      </div>
      <div className="settings-content">
        <Form
          form={form}
          layout="vertical"
        >
          <PipelineSettingsTabs />
        </Form>
      </div>
    </div>
  )
}
