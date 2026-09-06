set -euo pipefail
python - <<'PY'
import json, textwrap
from pathlib import Path
p=Path('/scratch/context_bundle.json')
data=json.loads(p.read_text())
for key in ['worker_brief','round','candidate','accepted_panel','current_panel','prompt','brief']:
    if key in data:
        print(f'--- {key} ---')
        print(json.dumps(data[key], indent=2)[:12000])
        print()
print('TOP KEYS:', list(data)[:50])
PY
echo '--- worker_brief.json if present ---'
if [ -f /scratch/worker_brief.json ]; then cat /scratch/worker_brief.json; else echo 'missing'; fi
