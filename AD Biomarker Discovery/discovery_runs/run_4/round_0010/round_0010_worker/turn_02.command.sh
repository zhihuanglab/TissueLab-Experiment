set -euo pipefail
python - <<'PY'
import json
from pprint import pprint
with open('/scratch/context_bundle.json') as f:
    ctx=json.load(f)
print("TOP-LEVEL KEYS:", list(ctx.keys()))
for k in ['worker_brief','round_plan','candidate','current_panel','accepted_panel','hypothesis','seed_candidate','round_context']:
    if k in ctx:
        print(f"\n--- {k} ---")
        pprint(ctx[k], width=120)
PY
