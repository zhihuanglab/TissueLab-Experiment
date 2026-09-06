set -euo pipefail
python - <<'PY'
from pathlib import Path
text=Path('/data/autoresearch_runs/run_20260421_001115_585b70/round_0006/round_0006_worker/result.py').read_text().splitlines()
for j in range(560, 582):
    print(text[j])
PY
