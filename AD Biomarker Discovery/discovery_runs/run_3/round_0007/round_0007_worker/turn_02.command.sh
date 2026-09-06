set -euo pipefail
python - <<'PY'
import json, os
with open('/scratch/context_bundle.json') as f:
    obj=json.load(f)
for key in ['worker_brief','accepted_panel','panel_state','round','candidate','hypothesis','brief']:
    if key in obj:
        print(f'--- {key} ---')
        print(json.dumps(obj[key], indent=2)[:12000])
        print()
print('top-level keys:', list(obj.keys()))
PY
echo '--- /scratch directory ---'
ls -1 /scratch | sed -n '1,120p'
if [ -f /scratch/worker_brief.json ]; then
  echo '--- worker_brief.json ---'
  python - <<'PY'
import json
with open('/scratch/worker_brief.json') as f:
    obj=json.load(f)
print(json.dumps(obj, indent=2)[:20000])
PY
fi
