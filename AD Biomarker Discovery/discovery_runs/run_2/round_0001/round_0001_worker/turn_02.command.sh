set -euo pipefail
python - <<'PY'
import json, pprint
with open('/scratch/context_bundle.json') as f:
    ctx=json.load(f)
print(ctx.keys())
for key in ['worker_brief','round','candidate','accepted_panel','current_panel','plan','hypothesis','current_accepted_panel']:
    if key in ctx:
        print(f"\n--- {key} ---")
        pprint.pp(ctx[key], width=120)
PY

echo
echo '--- /scratch/worker_brief.json (if present) ---'
if [ -f /scratch/worker_brief.json ]; then cat /scratch/worker_brief.json; else echo 'missing'; fi
