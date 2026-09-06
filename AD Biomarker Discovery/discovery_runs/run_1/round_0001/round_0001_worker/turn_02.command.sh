set -euo pipefail
python - <<'PY'
import json
with open('/scratch/context_bundle.json') as f:
    data=json.load(f)
print('keys:', sorted(data.keys()))
for k in ['worker_brief','current_panel','round_plan','candidate','hypothesis','accepted_panel','prior_rounds','brief']:
    if k in data:
        print(f'\n--- {k} ---')
        import pprint; pprint.pp(data[k])
PY
