set -euo pipefail
python - <<'PY'
from pathlib import Path
p=Path('/data/autoresearch_runs/run_20260421_001115_585b70/round_0009/round_0009_worker/result.py')
txt=p.read_text()
print(txt[16000:22000])
PY
