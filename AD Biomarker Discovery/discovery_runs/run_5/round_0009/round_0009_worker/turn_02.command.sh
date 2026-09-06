set -euo pipefail
python - <<'PY'
import json
from pathlib import Path
p=Path('/scratch/context_bundle.json')
obj=json.loads(p.read_text())
print('keys:', list(obj.keys()))
for k in ['worker_brief','round','candidate','accepted_panel','panel','hypothesis','current_panel','brief']:
    if k in obj:
        print(f'\n--- {k} ---')
        import pprint; pprint.pp(obj[k], width=120)
PY
echo
echo '--- worker_brief.json if present ---'
if [ -f /scratch/worker_brief.json ]; then cat /scratch/worker_brief.json; fi
