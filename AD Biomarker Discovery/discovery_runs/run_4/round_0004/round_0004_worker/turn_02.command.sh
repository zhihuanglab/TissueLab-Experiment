set -euo pipefail
python - <<'PY'
import json
from pathlib import Path
obj=json.loads(Path('/scratch/context_bundle.json').read_text())
print('keys:', sorted(obj))
for k in ['worker_brief','candidate','round','accepted_panel','current_panel','plan','hypothesis','task','brief']:
    if k in obj:
        print(f'\n--- {k} ---')
        import pprint; pprint.pp(obj[k])
PY
