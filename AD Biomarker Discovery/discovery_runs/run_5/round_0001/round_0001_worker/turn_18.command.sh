set -euo pipefail
python - <<'PY'
import json
payload=json.load(open('/scratch/results.json'))
for r in payload['ranked_variations']:
    print(r)
PY
