set -euo pipefail
python - <<'PY'
import json, pprint
obj=json.load(open('/scratch/context_bundle.json'))
rm=obj['runtime_manifest']
print(rm.keys())
print('slide manifest head', obj['slide_manifest_head'][:3] if isinstance(obj['slide_manifest_head'], list) else obj['slide_manifest_head'])
PY
