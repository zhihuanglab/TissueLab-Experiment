set -euo pipefail
python - <<'PY'
import json
from pathlib import Path
p=Path('/scratch/context_bundle.json')
obj=json.loads(p.read_text())
print('keys:', list(obj.keys()))
for k in obj:
    if isinstance(obj[k], dict):
        print('\n##', k, 'dict keys:', list(obj[k].keys())[:50])
    elif isinstance(obj[k], list):
        print('\n##', k, 'list len', len(obj[k]))
    else:
        print('\n##', k, type(obj[k]).__name__, obj[k] if type(obj[k]) in (str,int,float,bool) else '')
print('\nworker_brief maybe:')
for key in ['worker_brief','round_brief','candidate','hypothesis','plan','accepted_panel','current_panel']:
    if key in obj:
        print('\nKEY', key)
        print(json.dumps(obj[key], indent=2)[:4000])
PY
if [ -f /scratch/worker_brief.json ]; then
  echo '--- worker_brief.json ---'
  cat /scratch/worker_brief.json
fi
