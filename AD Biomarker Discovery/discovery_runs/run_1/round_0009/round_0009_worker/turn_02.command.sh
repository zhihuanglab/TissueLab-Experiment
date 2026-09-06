python - <<'PY'
import json
from pathlib import Path
obj=json.loads(Path('/scratch/context_bundle.json').read_text())
print(obj.keys())
for k in ['worker_brief','round_hypothesis','candidate','current_panel','accepted_panel','plan','round_plan','hypothesis_family','brief']:
    if k in obj:
        print(f'\n--- {k} ---')
        import pprint; pprint.pp(obj[k])
PY
