import { Dropdown, Button } from 'antd'
import type { MenuProps } from 'antd'
import { EllipsisOutlined } from '@ant-design/icons'

interface ActionItem {
  key: string
  label: string
  danger?: boolean
  disabled?: boolean
  onClick: () => void
}

interface TreeItemActionsProps {
  items: ActionItem[]
}

export default function TreeItemActions({ items }: TreeItemActionsProps) {
  const menuItems: NonNullable<MenuProps['items']> = items.map((item) => ({
    key: item.key,
    label: item.label,
    danger: item.danger,
    disabled: item.disabled,
  }))

  const onMenuClick: MenuProps['onClick'] = ({ key, domEvent }) => {
    domEvent.stopPropagation()
    const target = items.find((item) => item.key === key)
    if (!target || target.disabled) return
    target.onClick()
  }

  return (
    <div
      className="tree-item-actions"
      onClick={(e) => e.stopPropagation()}
      onMouseDown={(e) => e.stopPropagation()}
    >
      <Dropdown
        trigger={['click']}
        menu={{ items: menuItems, onClick: onMenuClick }}
      >
        <Button
          type="text"
          size="small"
          className="tree-action-btn"
          icon={<EllipsisOutlined />}
        />
      </Dropdown>
    </div>
  )
}
