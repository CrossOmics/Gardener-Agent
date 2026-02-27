import net from 'net'
import path from 'path'
import { app } from 'electron'
import { spawn, type ChildProcess } from 'child_process'
import fs from 'fs'

interface ServiceSpec {
  name: string
  port: number
  readyPath: string
  devExePathWin: string
  devExePathMac: string
  devAppExecutableMac: string
  devAppExecutableMacAlt: string
  packagedRelativeExePathWin: string
  packagedRelativeExePathMac: string
  packagedRelativeAppExecutableMac: string
  packagedRelativeAppExecutableMacAlt: string
  startupTimeoutMs: number
}

const SERVICE_SPECS: ServiceSpec[] = [
  {
    name: 'backend',
    port: 41888,
    readyPath: '/api/v1/project/root',
    devExePathWin: path.resolve(process.cwd(), 'release', 'lotus-backend.exe'),
    devExePathMac: path.resolve(process.cwd(), 'release', 'lotus-backend'),
    devAppExecutableMac: path.resolve(
      process.cwd(),
      'release',
      'lotus-backend.app',
      'Contents',
      'MacOS',
      'lotus-backend'
    ),
    devAppExecutableMacAlt: path.resolve(
      process.cwd(),
      'release',
      'LotusBackend.app',
      'Contents',
      'MacOS',
      'lotus-backend'
    ),
    packagedRelativeExePathWin: path.join('services', 'backend', 'lotus-backend.exe'),
    packagedRelativeExePathMac: path.join('services', 'backend', 'lotus-backend'),
    packagedRelativeAppExecutableMac: path.join(
      'services',
      'backend',
      'lotus-backend.app',
      'Contents',
      'MacOS',
      'lotus-backend'
    ),
    packagedRelativeAppExecutableMacAlt: path.join(
      'services',
      'backend',
      'LotusBackend.app',
      'Contents',
      'MacOS',
      'lotus-backend'
    ),
    startupTimeoutMs: 600000,
  },
  {
    name: 'agent',
    port: 41889,
    readyPath: '/docs',
    devExePathWin: path.resolve(process.cwd(), 'release', 'lotus-agent.exe'),
    devExePathMac: path.resolve(process.cwd(), 'release', 'lotus-agent'),
    devAppExecutableMac: path.resolve(
      process.cwd(),
      'release',
      'lotus-agent.app',
      'Contents',
      'MacOS',
      'lotus-agent'
    ),
    devAppExecutableMacAlt: path.resolve(
      process.cwd(),
      'release',
      'LotusAgent.app',
      'Contents',
      'MacOS',
      'lotus-agent'
    ),
    packagedRelativeExePathWin: path.join('services', 'agent', 'lotus-agent.exe'),
    packagedRelativeExePathMac: path.join('services', 'agent', 'lotus-agent'),
    packagedRelativeAppExecutableMac: path.join(
      'services',
      'agent',
      'lotus-agent.app',
      'Contents',
      'MacOS',
      'lotus-agent'
    ),
    packagedRelativeAppExecutableMacAlt: path.join(
      'services',
      'agent',
      'LotusAgent.app',
      'Contents',
      'MacOS',
      'lotus-agent'
    ),
    startupTimeoutMs: 600000,
  },
]

interface RunningService {
  spec: ServiceSpec
  child: ChildProcess
}

function isPortOpen(port: number): Promise<boolean> {
  return new Promise((resolve) => {
    const socket = new net.Socket()
    socket.setTimeout(400)

    socket.once('connect', () => {
      socket.destroy()
      resolve(true)
    })
    socket.once('timeout', () => {
      socket.destroy()
      resolve(false)
    })
    socket.once('error', () => {
      socket.destroy()
      resolve(false)
    })

    socket.connect(port, '127.0.0.1')
  })
}

async function waitForPort(port: number, timeoutMs: number): Promise<void> {
  const start = Date.now()
  while (Date.now() - start < timeoutMs) {
    if (await isPortOpen(port)) return
    await new Promise((r) => setTimeout(r, 300))
  }
  throw new Error(`Timed out waiting for service port ${port} after ${timeoutMs}ms`)
}

async function isHttpReady(spec: ServiceSpec): Promise<boolean> {
  try {
    const response = await fetch(`http://127.0.0.1:${spec.port}${spec.readyPath}`, {
      method: 'GET',
    })
    return response.ok
  } catch {
    return false
  }
}

async function waitForHttpReady(spec: ServiceSpec, timeoutMs: number): Promise<void> {
  const start = Date.now()
  while (Date.now() - start < timeoutMs) {
    if (await isHttpReady(spec)) return
    await new Promise((r) => setTimeout(r, 400))
  }
  throw new Error(
    `Timed out waiting for ${spec.name} HTTP readiness at http://127.0.0.1:${spec.port}${spec.readyPath} after ${timeoutMs}ms`
  )
}

function waitForChildExit(child: ChildProcess): Promise<never> {
  return new Promise((_, reject) => {
    child.once('error', (err) => {
      reject(new Error(`Service process failed to start: ${err.message}`))
    })
    child.once('exit', (code, signal) => {
      reject(new Error(`Service process exited before ready (code=${code ?? 'null'}, signal=${signal ?? 'null'})`))
    })
  })
}

function pickExistingPath(candidates: string[]): string | null {
  for (const candidate of candidates) {
    if (fs.existsSync(candidate)) return candidate
  }
  return null
}

function resolveExePath(spec: ServiceSpec): string {
  if (process.platform === 'win32') {
    const candidates = app.isPackaged
      ? [path.join(process.resourcesPath, spec.packagedRelativeExePathWin)]
      : [spec.devExePathWin]
    return pickExistingPath(candidates) ?? candidates[0]
  }

  if (process.platform === 'darwin') {
    const candidates = app.isPackaged
      ? [
          path.join(process.resourcesPath, spec.packagedRelativeExePathMac),
          path.join(process.resourcesPath, spec.packagedRelativeAppExecutableMac),
          path.join(process.resourcesPath, spec.packagedRelativeAppExecutableMacAlt),
        ]
      : [spec.devExePathMac, spec.devAppExecutableMac, spec.devAppExecutableMacAlt]
    return pickExistingPath(candidates) ?? candidates[0]
  }

  throw new Error(`Unsupported platform for local services: ${process.platform}`)
}

function stopProcessTree(child: ChildProcess): Promise<void> {
  return new Promise((resolve) => {
    if (!child.pid) {
      resolve()
      return
    }

    if (process.platform === 'win32') {
      const killer = spawn('taskkill', ['/pid', String(child.pid), '/t', '/f'], { windowsHide: true })
      killer.once('close', () => resolve())
      killer.once('error', () => resolve())
      return
    }

    try {
      child.kill('SIGTERM')
    } catch {
      // ignore
    }
    resolve()
  })
}

export function createLocalServiceManager() {
  const running = new Map<string, RunningService>()

  const startService = async (serviceName: string) => {
    const spec = SERVICE_SPECS.find((item) => item.name === serviceName)
    if (!spec) {
      throw new Error(`Unknown service: ${serviceName}`)
    }

    const existing = running.get(spec.name)
    if (existing) {
      await waitForHttpReady(spec, spec.startupTimeoutMs)
      return
    }

    const alreadyUp = await isPortOpen(spec.port)
    if (alreadyUp) {
      await waitForHttpReady(spec, spec.startupTimeoutMs)
      console.log(`[local-service] ${spec.name} already running and ready on :${spec.port}`)
      return
    }

    const exePath = resolveExePath(spec)
    if (!fs.existsSync(exePath)) {
      throw new Error(`Service executable not found: ${exePath}`)
    }
    console.log(`[local-service] starting ${spec.name} from ${exePath}`)
    const child = spawn(exePath, [], {
      cwd: path.dirname(exePath),
      windowsHide: true,
    })

    child.stdout?.on('data', (buf) => {
      console.log(`[${spec.name}] ${String(buf).trimEnd()}`)
    })
    child.stderr?.on('data', (buf) => {
      console.error(`[${spec.name}] ${String(buf).trimEnd()}`)
    })
    child.once('exit', (code, signal) => {
      console.log(`[local-service] ${spec.name} exited code=${code ?? 'null'} signal=${signal ?? 'null'}`)
      running.delete(spec.name)
    })

    running.set(spec.name, { spec, child })
    await Promise.race([
      waitForPort(spec.port, spec.startupTimeoutMs),
      waitForChildExit(child),
    ])
    await waitForHttpReady(spec, spec.startupTimeoutMs)
    console.log(`[local-service] ${spec.name} ready on :${spec.port}`)
  }

  const startAll = async () => {
    for (const spec of SERVICE_SPECS) {
      await startService(spec.name)
    }
  }

  const stopAll = async () => {
    const list = [...running.values()]
    running.clear()
    for (const item of list) {
      await stopProcessTree(item.child)
    }
  }

  return {
    startService,
    startAll,
    stopAll,
  }
}
