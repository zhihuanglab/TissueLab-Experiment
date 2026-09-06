set -euo pipefail
python - <<'PY'
import json
from pathlib import Path
obj=json.loads(Path('/scratch/context_bundle.json').read_text())
print('keys:', sorted(obj))
for k in sorted(obj):
    v=obj[k]
    if isinstance(v,(str,int,float,bool)) or v is None:
        print(f'\n## {k}\n{v}')
    elif isinstance(v,list):
        print(f'\n## {k} (list len={len(v)})')
        for i,item in enumerate(v[:5]):
            print(f'[{i}] {item}')
        if len(v)>5: print('...')
    elif isinstance(v,dict):
        print(f'\n## {k} (dict keys={list(v)[:20]})')
        # print selected subkeys fully if interesting
        if k in ('worker_brief','round','accepted_panel','candidate','plan'):
            import pprint; pprint.pp(v)
PY
