python - <<'PY'
from pathlib import Path
path=Path('/scratch/result.py')
text=path.read_text()
text=text.replace('"\\\\n".join(other_lines)', '"\\n".join(other_lines)')
text=text.replace('"\\\\n".join(lines)', '"\\n".join(lines)')
path.write_text(text)
print('reduced escapes')
PY
nl -ba /scratch/result.py | sed -n '484,558p'
