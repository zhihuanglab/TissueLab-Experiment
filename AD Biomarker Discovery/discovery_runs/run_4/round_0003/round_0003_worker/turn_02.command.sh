set -e
python - <<'PY'
import json, pprint
with open('/scratch/context_bundle.json') as f:
    obj=json.load(f)
print(obj.keys())
for k,v in obj.items():
    if k not in ['dataset_guide_excerpt','runtime_quickstart_excerpt']:
        print('\n###',k)
        pprint.pp(v, width=120)
PY
