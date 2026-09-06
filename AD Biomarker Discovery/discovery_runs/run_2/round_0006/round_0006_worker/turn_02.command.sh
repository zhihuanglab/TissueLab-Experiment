set -euo pipefail
python - <<'PY'
import json
from pprint import pprint
with open('/scratch/context_bundle.json') as f:
    obj=json.load(f)
print('TOP KEYS:', list(obj.keys()))
for k in ['worker_brief','current_panel','accepted_panel','round','candidate','hypothesis','plan','brief']:
    if k in obj:
        print(f'\n## {k}')
        pprint(obj[k])
print('\nLooking for strings mentioning candidate/family/variation...')
def walk(o,path=''):
    if isinstance(o,dict):
        for k,v in o.items():
            p=f'{path}.{k}' if path else k
            if isinstance(v,(dict,list)):
                walk(v,p)
            else:
                s=str(v)
                if any(tok in s.lower() for tok in ['variation','family','candidate','hypothesis','round','population','region','niche']):
                    print(p,':',s)
    elif isinstance(o,list):
        for i,v in enumerate(o):
            walk(v,f'{path}[{i}]')
walk(obj)
PY
[ -f /scratch/worker_brief.json ] && { echo '--- worker_brief.json ---'; sed -n '1,260p' /scratch/worker_brief.json; } || true
