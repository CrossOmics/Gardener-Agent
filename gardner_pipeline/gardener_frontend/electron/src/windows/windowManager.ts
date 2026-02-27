import { BrowserWindow } from 'electron'
import type { ProjectInfo } from '../types'
import { pathToFileURL } from 'url'
import fs from 'fs'

interface WindowManagerDeps {
  preloadPath: string
  getProjects: () => ProjectInfo[]
  rendererDevUrl: string
  rendererProdIndexPath: string
  isPackaged: boolean
}

export function createWindowManager({
  preloadPath,
  getProjects,
  rendererDevUrl,
  rendererProdIndexPath,
  isPackaged,
}: WindowManagerDeps) {
  const projectWindows = new Map<string, BrowserWindow>()
  const minWindowWidth = 1750
  const minWindowHeight = 1050

  const loadRenderer = (win: BrowserWindow, query: Record<string, string>) => {
    if (!isPackaged) {
      const url = new URL(rendererDevUrl)
      for (const [k, v] of Object.entries(query)) {
        url.searchParams.set(k, v)
      }
      void win.loadURL(url.toString())
      return
    }
    const fileUrl = pathToFileURL(rendererProdIndexPath)
    for (const [k, v] of Object.entries(query)) {
      fileUrl.searchParams.set(k, v)
    }
    if (!fs.existsSync(rendererProdIndexPath)) {
      console.error(`[window] renderer index not found: ${rendererProdIndexPath}`)
    }
    void win.loadURL(fileUrl.toString())
  }

  const createProjectWindow = (projectId: string) => {
    const existingWindow = projectWindows.get(projectId)
    if (existingWindow && !existingWindow.isDestroyed()) {
      existingWindow.focus()
      return existingWindow
    }

    const project = getProjects().find(p => p.id === projectId)
    const projectTitle = project?.name || projectId

    const win = new BrowserWindow({
      width: minWindowWidth,
      height: minWindowHeight,
      minWidth: minWindowWidth,
      minHeight: minWindowHeight,
      webPreferences: {
        preload: preloadPath,
      },
      title: `Lotus - ${projectTitle}`,
    })

    win.webContents.on('did-fail-load', (_event, code, desc, url) => {
      console.error(`[window] did-fail-load code=${code} desc=${desc} url=${url}`)
    })

    loadRenderer(win, { projectId })

    projectWindows.set(projectId, win)
    win.on('closed', () => {
      projectWindows.delete(projectId)
    })

    return win
  }

  const openProjectFromCurrentContext = (projectId: string) => {
    const focusedWindow = BrowserWindow.getFocusedWindow()
    const project = getProjects().find(p => p.id === projectId)
    const projectTitle = project?.name || projectId

    if (focusedWindow && !focusedWindow.isDestroyed()) {
      const currentUrl = focusedWindow.webContents.getURL()
      if (currentUrl.includes('welcome=true')) {
        focusedWindow.setTitle(`Lotus - ${projectTitle}`)
        loadRenderer(focusedWindow, { projectId })
        return focusedWindow
      }
    }

    return createProjectWindow(projectId)
  }

  const createWelcomeWindow = () => {
    const win = new BrowserWindow({
      width: minWindowWidth,
      height: minWindowHeight,
      minWidth: minWindowWidth,
      minHeight: minWindowHeight,
      webPreferences: {
        preload: preloadPath,
      },
      title: 'Lotus - Welcome',
    })

    win.webContents.on('did-fail-load', (_event, code, desc, url) => {
      console.error(`[window] did-fail-load code=${code} desc=${desc} url=${url}`)
    })

    loadRenderer(win, { welcome: 'true' })
    return win
  }

  return {
    createProjectWindow,
    openProjectFromCurrentContext,
    createWelcomeWindow,
  }
}
