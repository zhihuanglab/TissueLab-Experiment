set -euo pipefail
python - <<'PY'
import json, pprint
with open('/scratch/context_bundle.json') as f:
    data=json.load(f)
pprint.pp(data['worker_brief'])
PY
