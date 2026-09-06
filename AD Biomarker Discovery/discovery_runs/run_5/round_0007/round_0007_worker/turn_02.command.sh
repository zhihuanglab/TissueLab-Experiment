set -euo pipefail
python - <<'PY'
import json
obj=json.load(open('/scratch/context_bundle.json'))
print('TOP KEYS:', list(obj.keys()))
for k in ['worker_brief','current_panel','accepted_panel','round','candidate','hypothesis','family','plan','context','panel_state']:
    if k in obj:
        print(f'\n## {k}')
        print(json.dumps(obj[k], indent=2)[:6000])
PY
