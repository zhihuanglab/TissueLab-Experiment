set -e
python - <<'PY'
import json, pprint, pathlib
path = pathlib.Path('/scratch/context_bundle.json')
data = json.loads(path.read_text())
print('TOP KEYS:', sorted(data.keys()))
for k in ['worker_brief','round_plan','accepted_panel','current_panel','candidate','hypothesis','brief']:
    if k in data:
        print(f'\n== {k} ==')
        pprint.pp(data[k], sort_dicts=False)
print('\nSearch for likely plan-bearing keys:')
for k,v in data.items():
    if isinstance(v, dict):
        if any(tok in k.lower() for tok in ['brief','plan','candidate','round','task','worker']):
            print(k, 'keys=', list(v.keys())[:20])
PY
if [ -f /scratch/worker_brief.json ]; then echo '--- worker_brief.json ---'; cat /scratch/worker_brief.json; fi
