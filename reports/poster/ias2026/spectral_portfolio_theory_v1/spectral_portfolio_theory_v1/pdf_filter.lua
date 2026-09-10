function Code(el)
  if FORMAT:match('latex') then
    return pandoc.RawInline('latex', '\\nolinkurl{' .. el.text .. '}')
  end
end
