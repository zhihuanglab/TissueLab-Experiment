set -euo pipefail
python - <<'PY'
import json, pprint
obj=json.load(open('/scratch/context_bundle.json'))
print('TOP-LEVEL KEYS:', list(obj))
for k in ['worker_brief','round','candidate','accepted_panel','current_panel','hypothesis','analysis_plan','brief']:
    if k in obj:
        print(f'\n=== {k} ===')
        pprint.pp(obj[k])
PY
echo '--- worker_brief.json ---'
python - <<'PY'
import json, pprint
obj=json.load(open('/scratch/worker_brief.json'))
pprint.pp(obj)
PY
echo '--- preseeded result.py ---'
sed -n '1,260p' /scratch/result.py
