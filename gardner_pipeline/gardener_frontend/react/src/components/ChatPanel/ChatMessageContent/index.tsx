import React from 'react'
import './styles.css'

interface ChatMessageContentProps {
  content: string
}

type InlineToken =
  | { type: 'text'; value: string }
  | { type: 'code'; value: string }
  | { type: 'bold'; value: string }
  | { type: 'italic'; value: string }
  | { type: 'link'; value: string; href: string }

function parseInline(text: string): InlineToken[] {
  const tokens: InlineToken[] = []
  let i = 0

  const pushText = (value: string) => {
    if (!value) return
    tokens.push({ type: 'text', value })
  }

  while (i < text.length) {
    const remaining = text.slice(i)

    const codeMatch = remaining.match(/^`([^`]+)`/)
    if (codeMatch) {
      tokens.push({ type: 'code', value: codeMatch[1] })
      i += codeMatch[0].length
      continue
    }

    const boldMatch = remaining.match(/^\*\*([^*]+)\*\*/)
    if (boldMatch) {
      tokens.push({ type: 'bold', value: boldMatch[1] })
      i += boldMatch[0].length
      continue
    }

    const italicMatch = remaining.match(/^\*([^*]+)\*/)
    if (italicMatch) {
      tokens.push({ type: 'italic', value: italicMatch[1] })
      i += italicMatch[0].length
      continue
    }

    const linkMatch = remaining.match(/^\[([^\]]+)\]\((https?:\/\/[^\s)]+)\)/)
    if (linkMatch) {
      tokens.push({ type: 'link', value: linkMatch[1], href: linkMatch[2] })
      i += linkMatch[0].length
      continue
    }

    let nextSpecial = remaining.length
    const indices = [
      remaining.indexOf('`'),
      remaining.indexOf('**'),
      remaining.indexOf('*'),
      remaining.indexOf('['),
    ].filter((n) => n >= 0)
    if (indices.length > 0) {
      nextSpecial = Math.min(...indices)
    }

    pushText(remaining.slice(0, nextSpecial || 1))
    i += nextSpecial || 1
  }

  return tokens
}

function renderInline(text: string, keyPrefix: string) {
  return parseInline(text).map((token, idx) => {
    const key = `${keyPrefix}-${idx}`
    if (token.type === 'code') return <code key={key}>{token.value}</code>
    if (token.type === 'bold') return <strong key={key}>{token.value}</strong>
    if (token.type === 'italic') return <em key={key}>{token.value}</em>
    if (token.type === 'link') {
      return (
        <a key={key} href={token.href} target="_blank" rel="noopener noreferrer">
          {token.value}
        </a>
      )
    }
    return <React.Fragment key={key}>{token.value}</React.Fragment>
  })
}

function renderBlocks(content: string) {
  const lines = content.replace(/\r\n/g, '\n').split('\n')
  const blocks: React.ReactNode[] = []
  let i = 0

  while (i < lines.length) {
    const line = lines[i]

    if (line.trim() === '') {
      i += 1
      continue
    }

    if (line.startsWith('```')) {
      const codeLines: string[] = []
      i += 1
      while (i < lines.length && !lines[i].startsWith('```')) {
        codeLines.push(lines[i])
        i += 1
      }
      if (i < lines.length && lines[i].startsWith('```')) i += 1
      blocks.push(
        <pre key={`pre-${i}`} className="chat-md-pre">
          <code>{codeLines.join('\n')}</code>
        </pre>
      )
      continue
    }

    const headingMatch = line.match(/^(#{1,6})\s+(.+)$/)
    if (headingMatch) {
      const level = Math.min(6, headingMatch[1].length)
      const text = headingMatch[2]
      const tagName = `h${level}`
      blocks.push(
        React.createElement(
          tagName,
          { key: `h-${i}`, className: `chat-md-h${level}` },
          renderInline(text, `h-${i}`)
        )
      )
      i += 1
      continue
    }

    if (/^>\s+/.test(line)) {
      const quoteLines: string[] = []
      while (i < lines.length && /^>\s+/.test(lines[i])) {
        quoteLines.push(lines[i].replace(/^>\s+/, ''))
        i += 1
      }
      blocks.push(
        <blockquote key={`q-${i}`} className="chat-md-quote">
          {quoteLines.map((q, idx) => <p key={`q-${i}-${idx}`}>{renderInline(q, `q-${i}-${idx}`)}</p>)}
        </blockquote>
      )
      continue
    }

    if (/^[-*]\s+/.test(line)) {
      const items: string[] = []
      while (i < lines.length && /^[-*]\s+/.test(lines[i])) {
        items.push(lines[i].replace(/^[-*]\s+/, ''))
        i += 1
      }
      blocks.push(
        <ul key={`ul-${i}`} className="chat-md-ul">
          {items.map((item, idx) => <li key={`ul-${i}-${idx}`}>{renderInline(item, `ul-${i}-${idx}`)}</li>)}
        </ul>
      )
      continue
    }

    if (/^\d+\.\s+/.test(line)) {
      const items: string[] = []
      while (i < lines.length && /^\d+\.\s+/.test(lines[i])) {
        items.push(lines[i].replace(/^\d+\.\s+/, ''))
        i += 1
      }
      blocks.push(
        <ol key={`ol-${i}`} className="chat-md-ol">
          {items.map((item, idx) => <li key={`ol-${i}-${idx}`}>{renderInline(item, `ol-${i}-${idx}`)}</li>)}
        </ol>
      )
      continue
    }

    const paragraph: string[] = [line]
    i += 1
    while (
      i < lines.length &&
      lines[i].trim() !== '' &&
      !lines[i].startsWith('```') &&
      !/^(#{1,6})\s+/.test(lines[i]) &&
      !/^>\s+/.test(lines[i]) &&
      !/^[-*]\s+/.test(lines[i]) &&
      !/^\d+\.\s+/.test(lines[i])
    ) {
      paragraph.push(lines[i])
      i += 1
    }

    blocks.push(
      <p key={`p-${i}`} className="chat-md-p">
        {renderInline(paragraph.join(' '), `p-${i}`)}
      </p>
    )
  }

  return blocks
}

export default function ChatMessageContent({ content }: ChatMessageContentProps) {
  return <div className="chat-md">{renderBlocks(content)}</div>
}
