set -euo pipefail
python - <<'PY'
import json
from pathlib import Path
obj=json.loads(Path('/scratch/context_bundle.json').read_text())
print('TOP-LEVEL KEYS:', list(obj.keys()))
for k in ['worker_brief','candidate','round','accepted_panel','current_panel','hypothesis','plan','seed_candidate']:
    if k in obj:
        print(f'\n--- {k} ---')
        import pprint; pprint.pp(obj[k])
PY
