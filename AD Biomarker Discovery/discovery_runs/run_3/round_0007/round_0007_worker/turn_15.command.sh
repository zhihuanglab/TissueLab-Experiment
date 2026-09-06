set -euo pipefail
python - <<'PY'
from pathlib import Path
text=Path('/data/autoresearch_runs/run_20260421_001115_585b70/round_0006/round_0006_worker/result.py').read_text().splitlines()
for i,l in enumerate(text):
    if 'ranked_variations' in l or 'sort(' in l:
        if 520 < i < 620:
            print(i+1,l)
PY
