import { useCallback, useEffect, useMemo } from 'react'
import type { NavItem } from '../components/ProjectSidebar'

export function useArrowNavigation(
  navList: NavItem[],
  currentNavItem: NavItem | null,
  onNavigate: (item: NavItem) => void
) {
  const currentIndex = useMemo(() => {
    if (!currentNavItem) return -1
    return navList.findIndex((item) => item.path === currentNavItem.path)
  }, [navList, currentNavItem])

  const hasPrev = currentIndex > 0
  const hasNext = currentIndex >= 0 && currentIndex < navList.length - 1

  const handlePrev = useCallback(() => {
    if (!hasPrev) return
    onNavigate(navList[currentIndex - 1])
  }, [hasPrev, onNavigate, navList, currentIndex])

  const handleNext = useCallback(() => {
    if (!hasNext) return
    onNavigate(navList[currentIndex + 1])
  }, [hasNext, onNavigate, navList, currentIndex])

  useEffect(() => {
    const handleKeyDown = (e: KeyboardEvent) => {
      if (e.target instanceof HTMLInputElement || e.target instanceof HTMLTextAreaElement) {
        return
      }

      if (e.key === 'ArrowLeft' && hasPrev) {
        e.preventDefault()
        handlePrev()
      } else if (e.key === 'ArrowRight' && hasNext) {
        e.preventDefault()
        handleNext()
      }
    }

    window.addEventListener('keydown', handleKeyDown)
    return () => window.removeEventListener('keydown', handleKeyDown)
  }, [hasPrev, hasNext, handlePrev, handleNext])

  return {
    currentIndex,
    hasPrev,
    hasNext,
    handlePrev,
    handleNext,
  }
}
