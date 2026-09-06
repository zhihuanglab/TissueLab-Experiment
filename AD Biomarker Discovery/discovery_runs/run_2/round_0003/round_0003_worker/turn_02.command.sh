set -euo pipefail
python - <<'PY'
import json, textwrap
with open('/scratch/context_bundle.json') as f:
    c=json.load(f)
print('TOP-LEVEL KEYS:', sorted(c.keys()))
for k in ['worker_brief','round','accepted_panel','current_panel','candidate','hypothesis','plan']:
    if k in c:
        print(f'\n--- {k} ---')
        print(json.dumps(c[k], indent=2)[:12000])
print('\n--- worker_brief.json ---')
with open('/scratch/worker_brief.json') as f:
    print(f.read())
print('\n--- result.py template ---')
with open('/scratch/result.py') as f:
    print(f.read())
PY
