set -euo pipefail
python - <<'PY'
from pathlib import Path
paths=list(Path('/data/autoresearch_runs').rglob('round_0008_worker/results.json'))
print('n', len(paths))
for p in paths[:3]:
    print(p)
    print(p.read_text()[:2000])
PY
