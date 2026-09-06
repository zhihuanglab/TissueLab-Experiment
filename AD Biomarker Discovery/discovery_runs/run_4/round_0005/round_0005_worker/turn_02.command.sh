set -euo pipefail
python - <<'PY'
import json
with open('/scratch/context_bundle.json') as f:
    obj=json.load(f)
print(obj.keys())
for k,v in obj.items():
    if k not in ('dataset_guide_excerpt','runtime_quickstart_excerpt'):
        print('\nKEY',k)
        if isinstance(v, (dict,list)):
            import pprint;pprint.pp(v,width=120,compact=False)
        else:
            print(v)
PY
