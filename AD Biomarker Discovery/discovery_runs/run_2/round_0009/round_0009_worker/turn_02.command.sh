set -euo pipefail
python - <<'PY'
import json, pprint
with open('/scratch/context_bundle.json') as f:
    obj=json.load(f)
print('KEYS:', list(obj))
for k in ['worker_brief','round_plan','accepted_panel','current_panel','candidate','hypothesis','brief']:
    if k in obj:
        print(f'\n--- {k} ---')
        pprint.pp(obj[k], width=120)
PY
