set -euo pipefail
python - <<'PY'
import json
from pprint import pprint
with open('/scratch/context_bundle.json') as f:
    obj=json.load(f)
print('TOP-LEVEL KEYS:', list(obj))
for k in ['worker_brief','candidate','round','accepted_panel','current_panel','panel','brief']:
    if k in obj:
        print(f'\n== {k} ==')
        pprint(obj[k])
print('\nSEARCH FOR hypothesis-like fields:')
def walk(o, path=''):
    if isinstance(o, dict):
        for k,v in o.items():
            p = f'{path}.{k}' if path else k
            if any(tok in k.lower() for tok in ['hyp','candidate','brief','plan','variation','panel','round']):
                print('\nPATH',p)
                if isinstance(v,(dict,list)):
                    from pprint import pformat
                    print(pformat(v)[:4000])
                else:
                    print(repr(v))
            walk(v,p)
    elif isinstance(o,list):
        for i,v in enumerate(o[:10]):
            walk(v,f'{path}[{i}]')
walk(obj)
PY
