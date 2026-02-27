import { useState, useCallback, useEffect, useRef } from 'react'
import { Form, Button, message, Typography, Steps } from 'antd'
import { PlayCircleOutlined, ArrowLeftOutlined, CheckCircleOutlined, LoadingOutlined } from '@ant-design/icons'
import {
  settingsToRequest,
  type PipelineSettings,
  fetchLatestPreferenceDefaults,
} from '../../services/pipelineConfig'
import { runPipeline } from '../../services/endpoints'
import PipelineSettingsTabs from './PipelineSettingsTabs'
import './styles.css'

const { Title, Text, Paragraph } = Typography

interface PipelineSetupProps {
  projectId: string
  datasetId: string
  datasetName: string
  onBack: () => void
  onComplete: (datasetId: string, datasetName: string) => void
}

type SetupStep = 'configure' | 'running'

function inferGroupByFromMethod(method: PipelineSettings['clustering_method'] | undefined): string {
  if (method === 'louvain' || method === 'cplearn' || method === 'leiden') return method
  return 'leiden'
}

export default function PipelineSetup({
  projectId,
  datasetId,
  datasetName,
  onBack,
  onComplete,
}: PipelineSetupProps) {
  const [form] = Form.useForm<PipelineSettings>()
  const [messageApi, contextHolder] = message.useMessage()
  const [step, setStep] = useState<SetupStep>('configure')
  const [currentPipelineStep, setCurrentPipelineStep] = useState(0)
  const [selectedAnnotationModels, setSelectedAnnotationModels] = useState<string[]>([])
  const clusteringMethod = Form.useWatch('clustering_method', form)
  const prevInferredGroupByRef = useRef<string>('leiden')

  useEffect(() => {
    let mounted = true
    void (async () => {
      try {
        const defaults = await fetchLatestPreferenceDefaults()
        if (!mounted) return
        form.setFieldsValue(defaults.settings)
        setSelectedAnnotationModels(defaults.selectedAnnotationModels)
        prevInferredGroupByRef.current = inferGroupByFromMethod(defaults.settings.clustering_method)
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

  const handleAnnotationSelectionChange = useCallback((selectedIds: string[]) => {
    setSelectedAnnotationModels((prev) => {
      if (prev.length === selectedIds.length && prev.every((item, idx) => item === selectedIds[idx])) {
        return prev
      }
      return selectedIds
    })
  }, [])

  const pipelineSteps = [
    'QC Filtering',
    'HVG Selection',
    'PCA Analysis',
    'Neighbors Graph',
    'Clustering',
    'DGE Analysis',
    'Annotation',
  ]

  const handleRunPipeline = useCallback(async () => {
    try {
      const values = await form.validateFields()

      // Merge annotation models into settings
      const fullSettings: PipelineSettings = {
        ...values,
        selected_annotation_models: selectedAnnotationModels,
      }

      setStep('running')
      setCurrentPipelineStep(0)

      // Simulate step progression
      const stepInterval = setInterval(() => {
        setCurrentPipelineStep(prev => {
          if (prev < pipelineSteps.length - 1) {
            return prev + 1
          }
          clearInterval(stepInterval)
          return prev
        })
      }, 400)

      // Make the API call
      const request = settingsToRequest(projectId, datasetId, fullSettings)
      const response = await runPipeline(request)

      clearInterval(stepInterval)
      setCurrentPipelineStep(pipelineSteps.length)

      console.log('Pipeline completed:', response)
      messageApi.success(response.msg || 'Pipeline completed successfully!')
      onComplete(datasetId, datasetName)
    } catch (error) {
      console.error('Pipeline failed:', error)
      const reason = error instanceof Error ? error.message : 'Unknown error'
      messageApi.error(`Pipeline execution failed: ${reason}`)
      setStep('configure')
    }
  }, [datasetId, datasetName, form, messageApi, onComplete, pipelineSteps, projectId, selectedAnnotationModels])

  const renderConfigureStep = () => {
    return (
      <div className="pipeline-setup-container">
        <div className="pipeline-setup-header">
          <Button
            type="text"
            icon={<ArrowLeftOutlined />}
            onClick={onBack}
            className="back-btn"
          >
            Back
          </Button>
          <div className="pipeline-setup-title">
            <Title level={3}>Configure Analysis Pipeline</Title>
            <Text type="secondary">Dataset: {datasetName}</Text>
          </div>
          <Button
            type="primary"
            icon={<PlayCircleOutlined />}
            size="large"
            onClick={handleRunPipeline}
            className="run-btn"
          >
            Run Pipeline
          </Button>
        </div>

        <div className="pipeline-setup-content">
          <Form
            form={form}
            layout="vertical"
          >
            <PipelineSettingsTabs
              showAnnotationModelSelector
              selectedAnnotationModels={selectedAnnotationModels}
              onAnnotationSelectionChange={handleAnnotationSelectionChange}
            />
          </Form>
        </div>
      </div>
    )
  }

  const renderRunningStep = () => (
    <div className="pipeline-running-container">
      <Title level={3}>Running Analysis Pipeline</Title>
      <Paragraph type="secondary">
        Processing dataset: {datasetName}
      </Paragraph>

      <div className="pipeline-steps">
        <Steps
          direction="vertical"
          current={currentPipelineStep}
          items={pipelineSteps.map((stepName, index) => ({
            title: stepName,
            icon: index < currentPipelineStep ? (
              <CheckCircleOutlined />
            ) : index === currentPipelineStep ? (
              <LoadingOutlined />
            ) : undefined,
          }))}
        />
      </div>

      <Paragraph type="secondary" style={{ marginTop: 24 }}>
        This may take a few minutes depending on your dataset size...
      </Paragraph>
    </div>
  )

  return (
    <div className="pipeline-setup-page">
      {contextHolder}
      {step === 'configure' && renderConfigureStep()}
      {step === 'running' && renderRunningStep()}
    </div>
  )
}
