set -euo pipefail
python - <<'PY'
import json
from pprint import pprint
with open('/scratch/context_bundle.json') as f:
    obj=json.load(f)
print('TOP LEVEL KEYS:', list(obj.keys()))
for k in ['worker_brief','accepted_panel','round','candidate','hypothesis','current_panel','plan']:
    if k in obj:
        print(f'\n## {k}')
        pprint(obj[k], sort_dicts=False)
PY
