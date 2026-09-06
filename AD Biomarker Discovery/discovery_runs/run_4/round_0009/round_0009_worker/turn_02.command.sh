set -euo pipefail
python - <<'PY'
import json
from pathlib import Path
p=Path('/scratch/context_bundle.json')
obj=json.loads(p.read_text())
print('Top-level keys:', list(obj.keys()))
for k in obj:
    if isinstance(obj[k], dict):
        print('\nKEY',k,'subkeys',list(obj[k].keys())[:50])
    elif isinstance(obj[k], list):
        print('\nKEY',k,'list len',len(obj[k]))
    else:
        print('\nKEY',k, type(obj[k]).__name__, obj[k])

print('\nPossible round-specific content:')
def find(node, path=''):
    if isinstance(node, dict):
        for kk,v in node.items():
            if any(tok in kk.lower() for tok in ['brief','candidate','hypothesis','round','panel','accepted']):
                print(path+'/'+kk, type(v).__name__)
                if not isinstance(v,(dict,list)):
                    print('  ',v)
            find(v, path+'/'+kk)
    elif isinstance(node, list):
        for i,v in enumerate(node[:20]):
            find(v, path+f'[{i}]')
find(obj)

for fname in ['/scratch/worker_brief.json','/scratch/worker_brief.md']:
    pp=Path(fname)
    if pp.exists():
        print(f'\n--- {fname} ---')
        print(pp.read_text()[:20000])
PY
