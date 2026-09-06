set -euo pipefail
python - <<'PY'
import json
obj=json.load(open('/scratch/context_bundle.json'))
print(obj.keys())
for k in ['worker_brief','round_hypothesis','candidate_family','accepted_panel','current_panel','brief','round_plan','hypothesis']:
    if k in obj:
        print(f"\n--- {k} ---")
        print(json.dumps(obj[k], indent=2)[:12000])
PY
