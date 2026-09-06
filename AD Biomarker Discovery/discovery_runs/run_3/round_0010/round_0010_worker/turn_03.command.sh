set -euo pipefail
python - <<'PY'
import json, pathlib
obj=json.loads(pathlib.Path('/scratch/context_bundle.json').read_text())
print(sorted(obj))
for k,v in obj.items():
    if isinstance(v, dict):
        print(k, 'dict keys=', sorted(v)[:20], '... total', len(v))
    elif isinstance(v, list):
        print(k, 'list len', len(v))
    else:
        print(k, type(v).__name__, (str(v)[:120]+'...') if len(str(v))>120 else str(v))
PY
