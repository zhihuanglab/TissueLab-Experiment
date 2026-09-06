set -euo pipefail
python - <<'PY'
import json
from pprint import pprint
with open('/scratch/context_bundle.json') as f:
    obj=json.load(f)
print('TOP KEYS:', list(obj.keys()))
for k in obj:
    if k not in ('dataset_guide_excerpt','runtime_quickstart_excerpt'):
        print(f'\n## {k}')
        pprint(obj[k], width=120)
PY
