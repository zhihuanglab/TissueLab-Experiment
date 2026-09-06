set -euo pipefail
python - <<'PY'
import json
from pprint import pprint
with open('/scratch/context_bundle.json') as f:
    data=json.load(f)
print('KEYS:', list(data.keys()))
for k in ['worker_brief','candidate','accepted_panel','round','hypothesis','plan','current_panel','panel']:
    if k in data:
        print(f'\n--- {k} ---')
        pprint(data[k], width=120)
PY
