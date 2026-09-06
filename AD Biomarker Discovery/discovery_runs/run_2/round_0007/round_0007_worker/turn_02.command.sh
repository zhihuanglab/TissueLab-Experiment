set -euo pipefail
python - <<'PY'
import json, pathlib, pprint
ctx=json.loads(pathlib.Path('/scratch/context_bundle.json').read_text())
print('TOP-LEVEL KEYS:', list(ctx.keys()))
for k in ['worker_brief','candidate','round','accepted_panel','current_panel','prior_rounds','notes']:
    if k in ctx:
        print(f'\n--- {k} ---')
        pprint.pp(ctx[k])
print('\n--- worker_brief.json ---')
wb=json.loads(pathlib.Path('/scratch/worker_brief.json').read_text())
pprint.pp(wb)
print('\n--- result.py template ---')
print(pathlib.Path('/scratch/result.py').read_text())
PY
