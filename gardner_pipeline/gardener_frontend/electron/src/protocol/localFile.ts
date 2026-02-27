import fs from 'fs'
import { protocol } from 'electron'

export function registerLocalFileSchemeAsPrivileged() {
  protocol.registerSchemesAsPrivileged([
    {
      scheme: 'local-file',
      privileges: {
        secure: true,
        standard: true,
        supportFetchAPI: true,
        stream: true,
      },
    },
  ])
}

export function registerLocalFileProtocol() {
  protocol.handle('local-file', (request) => {
    const url = new URL(request.url)
    const decodedPath = decodeURIComponent(url.pathname)
    let filePath = ''

    if (process.platform === 'win32') {
      if (url.hostname) {
        // local-file://C:/path/to/file.png -> C:\path\to\file.png
        filePath = `${url.hostname.toUpperCase()}:${decodedPath}`
      } else if (/^\/[a-zA-Z]:\//.test(decodedPath)) {
        // local-file:///C:/path/to/file.png -> C:\path\to\file.png
        filePath = decodedPath.slice(1)
      } else {
        filePath = decodedPath
      }
      filePath = filePath.replace(/\//g, '\\')
    } else {
      // macOS/Linux keep POSIX path. Also recover malformed URLs like:
      // local-file://users/name/file.png -> /users/name/file.png
      if (url.hostname) {
        filePath = `/${url.hostname}${decodedPath}`
      } else {
        filePath = decodedPath
      }
    }

    if (!fs.existsSync(filePath)) {
      console.error('File not found:', filePath)
      return new Response('File not found', { status: 404 })
    }

    const buffer = fs.readFileSync(filePath)
    const mimeType = filePath.toLowerCase().endsWith('.png') ? 'image/png' : 'application/octet-stream'

    return new Response(buffer, {
      headers: { 'Content-Type': mimeType },
    })
  })
}
