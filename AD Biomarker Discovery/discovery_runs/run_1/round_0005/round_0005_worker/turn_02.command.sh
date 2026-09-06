set -e
echo '--- worker_brief.json ---'
python - <<'PY'
import json, pprint
with open('/scratch/worker_brief.json') as f:
    obj=json.load(f)
pprint.pp(obj)
PY
echo '--- context keys and round info ---'
python - <<'PY'
import json
with open('/scratch/context_bundle.json') as f:
    obj=json.load(f)
print('keys:', list(obj.keys()))
for k in ['worker_brief','candidate','current_panel','accepted_panel','round','round_id','hypothesis','candidate_family','plan']:
    if k in obj:
        print(f'\n## {k}')
        print(obj[k] if not isinstance(obj[k], dict) else {kk: obj[k][kk] for kk in list(obj[k].keys())[:20]})
PY
echo '--- seeded result.py template ---'
sed -n '1,240p' /scratch/result.py
