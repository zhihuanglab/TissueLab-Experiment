set -euo pipefail
python - <<'PY'
import json
from pprint import pprint
with open('/scratch/context_bundle.json') as f:
    ctx=json.load(f)
print("TOP KEYS:", list(ctx.keys()))
for k in ['worker_brief','round_plan','accepted_panel','current_panel','candidate','round','brief','hypothesis','seed_candidate','panel_state']:
    if k in ctx:
        print(f"\n--- {k} ---")
        pprint(ctx[k], width=120)
PY
