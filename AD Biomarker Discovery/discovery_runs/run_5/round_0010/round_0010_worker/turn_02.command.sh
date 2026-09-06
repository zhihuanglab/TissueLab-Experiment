set -euo pipefail
python - <<'PY'
import json, textwrap
from pathlib import Path
p=Path('/scratch/context_bundle.json')
data=json.loads(p.read_text())
print('top keys:', list(data.keys()))
for k in ['worker_brief','candidate','round','accepted_panel','current_panel','hypothesis','plan','task']:
    if k in data:
        print(f'\n--- {k} ---')
        print(json.dumps(data[k], indent=2)[:4000])
PY
echo
echo '--- worker_brief.json if present ---'
if [ -f /scratch/worker_brief.json ]; then cat /scratch/worker_brief.json; else echo 'absent'; fi
