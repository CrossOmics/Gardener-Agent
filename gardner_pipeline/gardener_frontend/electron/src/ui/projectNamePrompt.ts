import { BrowserWindow, ipcMain } from 'electron'

export function promptForProjectName(): Promise<string | null> {
  return new Promise((resolve) => {
    const parentWindow = BrowserWindow.getFocusedWindow() || undefined
    const promptWindow = new BrowserWindow({
      width: 420,
      height: 220,
      parent: parentWindow,
      modal: !!parentWindow,
      resizable: false,
      minimizable: false,
      maximizable: false,
      autoHideMenuBar: true,
      title: 'Create New Project',
      webPreferences: {
        nodeIntegration: true,
        contextIsolation: false,
      },
    })

    let settled = false
    const finish = (value: string | null) => {
      if (settled) return
      settled = true
      ipcMain.removeListener('project-name-submitted', onSubmit)
      ipcMain.removeListener('project-name-cancelled', onCancel)
      resolve(value)
      if (!promptWindow.isDestroyed()) {
        promptWindow.close()
      }
    }

    const onSubmit = (_event: Electron.IpcMainEvent, value: string) => {
      const name = (value || '').trim()
      if (!name) {
        return
      }
      finish(name)
    }

    const onCancel = () => {
      finish(null)
    }

    ipcMain.once('project-name-submitted', onSubmit)
    ipcMain.once('project-name-cancelled', onCancel)

    promptWindow.on('closed', () => {
      finish(null)
    })

    const html = `
      <!doctype html>
      <html>
        <head>
          <meta charset="UTF-8" />
          <title>Create New Project</title>
          <style>
            body { font-family: -apple-system, Segoe UI, sans-serif; margin: 0; background: #ffffff; }
            .wrap { padding: 16px; }
            .title { font-size: 16px; font-weight: 600; margin-bottom: 12px; }
            .input { width: 100%; box-sizing: border-box; padding: 8px 10px; font-size: 14px; border: 1px solid #d9d9d9; border-radius: 6px; }
            .actions { display: flex; justify-content: flex-end; gap: 8px; margin-top: 14px; }
            .btn { padding: 6px 12px; font-size: 13px; border: 1px solid #d9d9d9; border-radius: 6px; background: #fff; cursor: pointer; }
            .btn.primary { background: #1677ff; color: #fff; border-color: #1677ff; }
          </style>
        </head>
        <body>
          <div class="wrap">
            <div class="title">Create New Project</div>
            <input id="project-name" class="input" type="text" placeholder="Enter project name" autofocus />
            <div class="actions">
              <button id="cancel-btn" class="btn" type="button">Cancel</button>
              <button id="create-btn" class="btn primary" type="button">Create</button>
            </div>
          </div>
          <script>
            const { ipcRenderer } = require('electron');
            const input = document.getElementById('project-name');
            const createBtn = document.getElementById('create-btn');
            const cancelBtn = document.getElementById('cancel-btn');

            const submit = () => {
              const value = (input.value || '').trim();
              if (!value) return;
              ipcRenderer.send('project-name-submitted', value);
            };

            createBtn.addEventListener('click', submit);
            cancelBtn.addEventListener('click', () => ipcRenderer.send('project-name-cancelled'));
            input.addEventListener('keydown', (e) => {
              if (e.key === 'Enter') submit();
              if (e.key === 'Escape') ipcRenderer.send('project-name-cancelled');
            });
          </script>
        </body>
      </html>
    `

    promptWindow.loadURL(`data:text/html;charset=utf-8,${encodeURIComponent(html)}`)
  })
}
