python - <<'PY'
import json, pathlib
print(pathlib.Path('/scratch/report.md').read_text().splitlines()[0:25])
data=json.load(open('/scratch/results.json'))
print(data['ci_lo'], data['ci_hi'])
PY
