import { Menu, dialog } from 'electron'
import type { ProjectInfo } from '../types'
import { promptForProjectName } from '../ui/projectNamePrompt'

interface BuildAppMenuDeps {
  getProjects: () => ProjectInfo[]
  createProjectViaApi: (projectName?: string) => Promise<ProjectInfo | null>
  openProjectFromCurrentContext: (projectId: string) => void
  createProjectWindow: (projectId: string) => void
}

export function buildAppMenu(deps: BuildAppMenuDeps): Menu {
  const projects = deps.getProjects()

  const template: Electron.MenuItemConstructorOptions[] = [
    {
      label: 'File',
      submenu: [
        {
          label: 'New Project',
          accelerator: 'CmdOrCtrl+N',
          click: async () => {
            const projectName = await promptForProjectName()
            if (!projectName) {
              return
            }

            const project = await deps.createProjectViaApi(projectName)
            if (project) {
              deps.openProjectFromCurrentContext(project.id)
              Menu.setApplicationMenu(buildAppMenu(deps))
            } else {
              dialog.showErrorBox('Create Project Failed', 'Failed to create project via backend API.')
            }
          },
        },
        { type: 'separator' },
        {
          label: 'Open Project',
          submenu: projects.length > 0
            ? projects.map(project => ({
                label: project.name,
                click: () => {
                  deps.createProjectWindow(project.id)
                },
              }))
            : [{ label: 'No projects available', enabled: false }],
        },
        { type: 'separator' },
        {
          label: 'Refresh Projects',
          accelerator: 'CmdOrCtrl+R',
          click: () => {
            Menu.setApplicationMenu(buildAppMenu(deps))
          },
        },
        { type: 'separator' },
        { role: 'quit' },
      ],
    },
    {
      label: 'Edit',
      submenu: [
        { role: 'undo' },
        { role: 'redo' },
        { type: 'separator' },
        { role: 'cut' },
        { role: 'copy' },
        { role: 'paste' },
      ],
    },
    {
      label: 'View',
      submenu: [
        { role: 'reload' },
        { role: 'forceReload' },
        { role: 'toggleDevTools' },
        { type: 'separator' },
        { role: 'resetZoom' },
        { role: 'zoomIn' },
        { role: 'zoomOut' },
        { type: 'separator' },
        { role: 'togglefullscreen' },
      ],
    },
    {
      label: 'Window',
      submenu: [
        { role: 'minimize' },
        { role: 'close' },
      ],
    },
  ]

  return Menu.buildFromTemplate(template)
}
