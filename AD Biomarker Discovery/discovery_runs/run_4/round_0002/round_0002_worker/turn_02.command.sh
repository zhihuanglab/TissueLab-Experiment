python - <<'PY'
import json
from pathlib import Path
p=Path('/scratch/context_bundle.json')
data=json.loads(p.read_text())
for k in data.keys():
    print('KEY',k, type(data[k]).__name__)
print('\n--- maybe brief-like keys ---')
for key in ['worker_brief','brief','candidate','round_plan','hypothesis','accepted_panel','current_panel','panel','task']:
    if key in data:
        print(f'\n## {key}')
        print(json.dumps(data[key], indent=2)[:20000])
PY
