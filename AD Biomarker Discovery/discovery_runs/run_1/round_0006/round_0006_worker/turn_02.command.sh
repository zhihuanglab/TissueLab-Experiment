set -euo pipefail
python - <<'PY'
import json
with open('/scratch/context_bundle.json') as f:
    ctx=json.load(f)
print('KEYS:', sorted(ctx.keys()))
for k in ['worker_brief','round_brief','candidate_brief','hypothesis','current_panel','accepted_panel']:
    if k in ctx:
        print(f'\n--- {k} ---')
        print(json.dumps(ctx[k], indent=2)[:12000])
PY

echo '--- worker_brief.json ---'
python - <<'PY'
import json
with open('/scratch/worker_brief.json') as f:
    obj=json.load(f)
print(json.dumps(obj, indent=2)[:20000])
PY
