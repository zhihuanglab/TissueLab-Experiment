set -euo pipefail
echo '--- worker_brief.json ---'
python - <<'PY'
import json, os
for p in ['/scratch/worker_brief.json','/scratch/context_bundle.json']:
    if os.path.exists(p):
        print(f'PATH {p}')
        with open(p) as f:
            obj=json.load(f)
        # print keys
        print('keys:', list(obj.keys())[:50])
        if p.endswith('worker_brief.json'):
            print(json.dumps(obj, indent=2)[:25000])
        else:
            # search recursively for strings containing hypothesis/candidate/round/family
            hits=[]
            def walk(x, path=''):
                if isinstance(x, dict):
                    for k,v in x.items():
                        walk(v, path + '/' + str(k))
                elif isinstance(x, list):
                    for i,v in enumerate(x):
                        walk(v, path + f'/{i}')
                elif isinstance(x, str):
                    s=x.lower()
                    if any(t in s for t in ['hypothesis','candidate','variation','family','round','accepted panel','biomarker']):
                        hits.append((path, x[:400]))
            walk(obj)
            for path, val in hits[:80]:
                print(path, '=>', val.replace('\n',' ')[:500])
PY
