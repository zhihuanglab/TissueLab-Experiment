set -euo pipefail
python - <<'PY'
import json, pathlib, pprint
obj=json.loads(pathlib.Path('/scratch/context_bundle.json').read_text())
print('KEYS:', sorted(obj.keys()))
for k in obj:
    if k not in ('dataset_guide_excerpt','runtime_quickstart_excerpt'):
        print(f'\n--- {k} ---')
        pprint.pp(obj[k], width=120)
PY
