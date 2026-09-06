set -euo pipefail
python - <<'PY'
from pathlib import Path
p=Path('/data/autoresearch_runs/run_20260421_001115_585b70/round_0008/round_0008_worker/results.json')
print(p.read_text()[:4000])
PY
