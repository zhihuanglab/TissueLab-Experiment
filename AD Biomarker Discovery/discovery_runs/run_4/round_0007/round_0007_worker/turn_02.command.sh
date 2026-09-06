set -euo pipefail
python - <<'PY'
import json, textwrap
with open('/scratch/context_bundle.json') as f:
    obj=json.load(f)
print("TOP-LEVEL KEYS:", sorted(obj.keys()))
for k in ['worker_brief','round','candidate','accepted_panel','panel','current_panel','hypothesis','brief']:
    if k in obj:
        print(f"\n--- {k} ---")
        print(json.dumps(obj[k], indent=2)[:12000])
PY
