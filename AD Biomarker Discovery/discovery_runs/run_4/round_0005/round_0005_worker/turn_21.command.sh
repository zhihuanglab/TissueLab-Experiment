set -euo pipefail
python - <<'PY'
from pathlib import Path
import json
print(Path('/scratch/result.py').exists(), Path('/scratch/report.md').exists(), Path('/scratch/results.json').exists())
res=json.load(open('/scratch/results.json'))
print(res['best_variation'], res['feature_column'])
print('result.py head:')
for i,line in enumerate(open('/scratch/result.py')):
    if i>25: break
    print(line.rstrip())
PY
