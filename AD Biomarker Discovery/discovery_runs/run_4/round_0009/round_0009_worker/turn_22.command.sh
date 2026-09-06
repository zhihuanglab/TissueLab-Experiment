python - <<'PY'
from pathlib import Path
path = Path('/scratch/result.py')
text = path.read_text()
text = text.replace('other_text = "\\n".join(other_lines) if other_lines else "- No alternative variations."', 'TEMP_PLACEHOLDER')  # if already fixed, unlikely
text = text.replace('other_text = "\n".join(other_lines) if other_lines else "- No alternative variations."', 'other_text = "\\\\n".join(other_lines) if other_lines else "- No alternative variations."')
text = text.replace('path.write_text("\n".join(lines), encoding="utf-8")', 'path.write_text("\\\\n".join(lines), encoding="utf-8")')
# cleanup if placeholder used
text = text.replace('TEMP_PLACEHOLDER', 'other_text = "\\\\n".join(other_lines) if other_lines else "- No alternative variations."')
path.write_text(text)
print('fixed literals')
PY
nl -ba /scratch/result.py | sed -n '482,560p'
