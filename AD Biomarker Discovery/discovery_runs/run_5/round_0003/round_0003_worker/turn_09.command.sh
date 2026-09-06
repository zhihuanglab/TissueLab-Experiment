set -euo pipefail
echo '--- round 1 result summary ---'
python - <<'PY'
import json, pathlib
for p in [
'/data/autoresearch_runs/run_20260421_052101_efbee1/round_0001/round_0001_worker/results.json',
'/data/autoresearch_runs/run_20260421_052101_efbee1/round_0002/round_0002_worker/results.json',
]:
    pp=pathlib.Path(p)
    if pp.exists():
        data=json.loads(pp.read_text())
        print(pp)
        for k in ['feature_name','feature_column','best_variation','partial_r','selection_score','loo_predictive_r','panel_candidate_score','delta_panel_score']:
            if k in data: print(' ',k,':',data[k])
        if 'ranked_variations' in data:
            print(' ranked head:', data['ranked_variations'][:3])
        print()
    else:
        print('missing', pp)
PY
