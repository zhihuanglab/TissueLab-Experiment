set -euo pipefail
python - <<'PY'
import json
from pathlib import Path
p=Path('/scratch/context_bundle.json')
data=json.loads(p.read_text())
print('KEYS:', list(data.keys()))
for k,v in data.items():
    if k in ('dataset_guide_excerpt','runtime_quickstart_excerpt'):
        continue
    print(f'\n## {k}')
    if isinstance(v, (dict,list)):
        import pprint; pprint.pp(v, width=120)
    else:
        print(v)
PY

echo '--- worker template ---'
sed -n '1,260p' /shared/templates/worker_analysis_template.py
