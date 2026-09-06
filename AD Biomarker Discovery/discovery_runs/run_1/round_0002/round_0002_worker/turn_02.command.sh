set -euo pipefail
python - <<'PY'
import json
from pathlib import Path
p=Path('/scratch/context_bundle.json')
obj=json.loads(p.read_text())
print("TOP KEYS:", list(obj))
for k in obj:
    if k not in ('dataset_guide_excerpt','runtime_quickstart_excerpt'):
        print(f"\n## {k}")
        v=obj[k]
        if isinstance(v,(dict,list)):
            import pprint; pprint.pp(v)
        else:
            print(v)
PY
