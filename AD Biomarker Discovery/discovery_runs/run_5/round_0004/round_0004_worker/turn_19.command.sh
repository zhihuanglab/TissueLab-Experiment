python - <<'PY'
from pathlib import Path
path=Path('/scratch/report.md')
text=path.read_text()
text=text.replace('bootstrap median partial r: **-0.3830**','bootstrap median partial r: **-0.3852**')
text=text.replace('bootstrap sign consistency: **0.9300**','bootstrap sign consistency: **0.9225**')
path.write_text(text)
PY
sed -n '1,220p' /scratch/report.md
