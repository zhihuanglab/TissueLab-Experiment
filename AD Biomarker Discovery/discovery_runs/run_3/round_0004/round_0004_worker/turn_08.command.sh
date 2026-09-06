set -euo pipefail
python - <<'PY'
from pathlib import Path
import json, os, glob
paths = sorted(Path('/data/autoresearch_runs').rglob('round_0003_worker/results.json'))
print(paths[:3])
if paths:
    p=paths[0]
    print('using',p)
    print(p.read_text()[:2000])
PY
