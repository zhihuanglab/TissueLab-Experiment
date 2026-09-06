set -euo pipefail
python - <<'PY'
import json, pprint
with open('/scratch/context_bundle.json') as f:
    obj=json.load(f)
print("keys:", list(obj.keys()))
for k in obj:
    if k not in ('dataset_guide_excerpt','runtime_quickstart_excerpt'):
        print("\nKEY",k)
        pprint.pp(obj[k], width=120)
PY
