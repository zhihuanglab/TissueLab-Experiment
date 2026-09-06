set -euo pipefail
python - <<'PY'
import json
with open('/scratch/context_bundle.json') as f:
    obj=json.load(f)
print('keys:', list(obj.keys()))
for k in ['worker_brief','round','candidate','accepted_panel','current_panel','hypothesis','plan','seed_candidate','prior_rounds_summary']:
    if k in obj:
        print(f'\n--- {k} ---')
        import pprint; pprint.pp(obj[k])
PY
