set -euo pipefail
python - <<'PY'
from pathlib import Path
text=Path('/data/autoresearch_runs/run_20260421_001115_585b70/round_0006/round_0006_worker/result.py').read_text().splitlines()
for i,l in enumerate(text):
    if 'ranked = sorted' in l:
        for j in range(i, i+12):
            print(text[j])
        break
PY
