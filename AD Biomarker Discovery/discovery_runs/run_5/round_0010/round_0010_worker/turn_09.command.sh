set -euo pipefail
python - <<'PY'
import json, glob, pathlib
for p in sorted(glob.glob('/Volumes/SSK SSD/TissueLab Revision/LFB/Seed48_Split/Training/autoresearch_runs/run_20260421_052101_efbee1/round_000*/round_*_worker/results.json'))[-3:]:
    data=json.load(open(p))
    print('\n---', p, '---')
    print('keys', list(data.keys())[:30])
    for k in ['best_variation','feature_name','feature_column','partial_r','selection_score','loo_predictive_r','is_loo_gap','adjusted_score']:
        if k in data:
            print(k, data[k])
PY
