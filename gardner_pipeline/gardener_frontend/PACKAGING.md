# Lotus Frontend Packaging Guide (ZH / EN)

## 1. Scope / 适用范围

This guide is for packaging the desktop app from:
- `lotus_frontend/react` (frontend static assets)
- `lotus_frontend/electron` (Electron shell)
- service binaries/apps placed under `lotus_frontend/electron/release`

本文档用于说明如何从以下目录打包桌面应用：
- `lotus_frontend/react`（前端静态资源）
- `lotus_frontend/electron`（Electron 壳层）
- 放在 `lotus_frontend/electron/release` 下的本地服务程序

---

## 2. Runtime Behavior / 运行逻辑

When the packaged app starts:
1. Electron starts backend and agent services.
2. It waits for ports `41888` (backend) and `41889` (agent).
3. It loads React assets from `resources/react/dist`.
4. On app quit, it stops both services.

打包后的应用启动时：
1. Electron 先启动 backend 和 agent 两个服务。
2. 等待端口 `41888`（backend）和 `41889`（agent）就绪。
3. 从 `resources/react/dist` 加载 React 页面。
4. 退出应用时自动停止两个服务。

---

## 3. Required Inputs / 必备输入

### Windows build inputs
- `electron/release/lotus-backend.exe`
- `electron/release/lotus-agent.exe`

### macOS build inputs
Current mac packaging config expects:
- `electron/release/LotusBackend.app`
- `electron/release/LotusAgent.app`

Each app must contain:
- `Contents/Info.plist`
- `Contents/MacOS/<executable>`

Executable paths should be:
- `LotusBackend.app/Contents/MacOS/lotus-backend`
- `LotusAgent.app/Contents/MacOS/lotus-agent`

Windows 打包输入：
- `electron/release/lotus-backend.exe`
- `electron/release/lotus-agent.exe`

macOS 打包输入（当前配置）：
- `electron/release/LotusBackend.app`
- `electron/release/LotusAgent.app`

---

## 4. Build Commands / 打包命令

Run commands in `lotus_frontend/electron`.

在 `lotus_frontend/electron` 目录执行命令。

### Windows installer (NSIS)
```powershell
npm run dist:win
```

### macOS DMG (default arch)
```bash
npm run dist:mac
```

### macOS DMG (Apple Silicon / M1/M2)
```bash
npm run dist:mac:arm64
```

These commands will:
- build React production assets
- compile Electron TypeScript
- run `electron-builder`

这些命令会：
- 构建 React 生产资源
- 编译 Electron TypeScript
- 调用 `electron-builder` 进行打包

---

## 5. Output Paths / 产物路径

### Windows
- Installer: `electron/release/Lotus Biological Agent Setup 1.0.0.exe`
- Unpacked app: `electron/release/win-unpacked/`

### macOS
- DMG output: `electron/release/*.dmg`

Windows 产物：
- 安装包：`electron/release/Lotus Biological Agent Setup 1.0.0.exe`
- 解包目录：`electron/release/win-unpacked/`

macOS 产物：
- DMG 文件：`electron/release/*.dmg`

---

## 6. Important Notes / 重要说明

1. `Setup.exe` / `.dmg` is the distributable artifact.
2. The app executable is not standalone; it depends on the bundled `resources/*`.
3. `.dmg` is an installer image, not a backend/agent binary.

1. `Setup.exe` / `.dmg` 是发布产物。
2. 主程序不是单文件独立运行，依赖 `resources/*`。
3. `.dmg` 是安装镜像，不是 backend/agent 可执行文件。

---

## 7. Troubleshooting / 常见问题

### A) PowerShell execution policy blocks npm scripts
```powershell
cmd /c npm run dist:win
```

### B) Startup failed: timed out waiting for service port
Check:
- service files exist under `electron/release`
- service can run manually
- ports are not occupied before app launch

### C) Access denied during packaging
Stop running app/services first:
```powershell
taskkill /IM "Lotus Biological Agent.exe" /F /T
taskkill /IM lotus-backend.exe /F /T
taskkill /IM lotus-agent.exe /F /T
```

### D) Blank window after launch
Ensure frontend build uses relative base (`base: './'` in Vite).

### E) How to run EXE with spaces in PowerShell
```powershell
& ".\Lotus Biological Agent.exe"
```

---

## 8. Release Checklist / 发布检查清单

1. `npm run dist:win` and/or `npm run dist:mac:arm64` succeeds.
2. Packaged app launches successfully.
3. Backend API works (`/api/v1/project/root` returns 200).
4. Agent API responds on port `41889`.
5. App quits cleanly and both service processes are terminated.

1. `npm run dist:win` 和/或 `npm run dist:mac:arm64` 成功。
2. 打包应用可正常启动。
3. Backend 接口可用（`/api/v1/project/root` 返回 200）。
4. Agent 在端口 `41889` 可响应。
5. 退出应用后服务进程能被正常关闭。
