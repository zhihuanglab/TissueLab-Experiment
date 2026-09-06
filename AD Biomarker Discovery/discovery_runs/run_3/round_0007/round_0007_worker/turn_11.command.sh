set -euo pipefail
python - <<'PY'
from pathlib import Path
for p in [
'/data/autoresearch_runs/run_20260421_001115_585b70/round_0006/round_0006_worker/result.py',
'/data/autoresearch_runs/run_20260421_001115_585b70/round_0005/round_0005_worker/result.py']:
    text=Path(p).read_text()
    for needle in ['adjusted_score','penalty']:
        idx=text.find(needle)
        print('FILE',p,'first',needle,'at',idx)
    print('---snip---')
    for line in text.splitlines():
        if 'adjusted_score' in line or 'penalty' in line:
            print(line)
    print('========')
PY
