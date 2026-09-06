set -euo pipefail
python - <<'PY'
import json
obj=json.load(open('/scratch/context_bundle.json'))
slides=obj['slide_manifest_head']['slides']
print(slides[0].keys())
for s in slides[:2]:
    print(s)
PY
