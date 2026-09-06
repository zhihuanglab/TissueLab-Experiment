set -euo pipefail
python - <<'PY'
import json, textwrap
with open('/scratch/context_bundle.json') as f:
    obj=json.load(f)
print('TOP KEYS:', sorted(obj.keys()))
for k in ['worker_brief','round_plan','candidate','current_panel','accepted_panel','worker_brief_path']:
    if k in obj:
        print(f'\n--- {k} ---')
        print(obj[k] if not isinstance(obj[k], (dict,list)) else json.dumps(obj[k], indent=2))
PY
