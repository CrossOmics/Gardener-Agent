import { Typography } from 'antd'
import { ProjectOutlined, DatabaseOutlined } from '@ant-design/icons'
import './styles.css'

const { Title, Paragraph } = Typography

interface ProjectSummaryProps {
  projectId: string
  projectName: string
}

export default function ProjectSummary({ projectId, projectName }: ProjectSummaryProps) {
  const summary = `Project "${projectName}" (${projectId}) is ready.`
  const datasetDescription = 'Import a dataset and run pipeline stages to populate analysis outputs.'

  return (
    <div className="project-summary-container">
      <div className="project-summary-header">
        <ProjectOutlined className="project-summary-icon" />
        <Title level={3} className="project-summary-title">{projectName}</Title>
      </div>

      <div className="project-summary-section">
        <Title level={5}>Summary</Title>
        <Paragraph className="project-summary-text">
          {summary}
        </Paragraph>
      </div>

      <div className="project-summary-section">
        <Title level={5}>
          <DatabaseOutlined style={{ marginRight: 8 }} />
          Dataset
        </Title>
        <Paragraph className="project-summary-text">
          {datasetDescription}
        </Paragraph>
      </div>
    </div>
  )
}
