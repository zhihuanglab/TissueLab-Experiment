set -euo pipefail
python - <<'PY'
import json
from pprint import pprint
with open('/scratch/context_bundle.json') as f:
    data=json.load(f)
print("TOP KEYS:", list(data.keys()))
for k in ['worker_brief','accepted_panel','round','candidate','hypothesis','plan','current_panel','panel_context','prior_rounds']:
    if k in data:
        print(f"\n--- {k} ---")
        pprint(data[k], width=120)
PY
