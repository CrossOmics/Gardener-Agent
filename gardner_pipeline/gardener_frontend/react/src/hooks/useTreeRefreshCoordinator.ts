import { useCallback, useRef } from 'react'
import type { MutableRefObject } from 'react'
import type { NavItem } from '../components/ProjectSidebar'

export type SidebarRefreshFn = (options?: { preferredNavItem?: NavItem | null }) => Promise<void>

interface UseTreeRefreshCoordinatorOptions {
  currentNavItem: NavItem | null
  sidebarRefreshRef: MutableRefObject<SidebarRefreshFn | null>
}

export function useTreeRefreshCoordinator({
  currentNavItem,
  sidebarRefreshRef,
}: UseTreeRefreshCoordinatorOptions) {
  const currentNavItemRef = useRef<NavItem | null>(currentNavItem)
  currentNavItemRef.current = currentNavItem

  const refreshTreeKeepingPosition = useCallback(async () => {
    const anchor = currentNavItemRef.current
    if (!sidebarRefreshRef.current) return
    await sidebarRefreshRef.current({ preferredNavItem: anchor })
  }, [sidebarRefreshRef])

  return {
    refreshTreeKeepingPosition,
  }
}
