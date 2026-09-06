set -euo pipefail
python - <<'PY'
import json
from pathlib import Path
p=Path('/scratch/context_bundle.json')
data=json.loads(p.read_text())
print('TOP KEYS:', list(data))
for k,v in data.items():
    if isinstance(v, (str,int,float,bool)) or v is None:
        print(f'\n## {k}\n{v if len(str(v))<3000 else str(v)[:3000]+"..."}')
    else:
        print(f'\n## {k} type={type(v).__name__}')
        if isinstance(v, dict):
            print('keys:', list(v)[:50])
            for sk,sv in v.items():
                if sk in ('worker_brief','candidate','round_plan','accepted_panel','current_panel','candidate_brief'):
                    print(f'-- {sk}: {sv}')
        elif isinstance(v, list):
            print('len=', len(v))
            for i,item in enumerate(v[:3]):
                print(i, item)
PY
