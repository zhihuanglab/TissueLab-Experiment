set -euo pipefail
echo '--- round 7 result.py ---'
sed -n '1,260p' /Volumes/SSK\ SSD/TissueLab\ Revision/LFB/Seed48_Split/Training/autoresearch_runs/run_20260420_215745_90b105/round_0007/round_0007_worker/result.py
echo
echo '--- round 7 results.json ---'
python - <<'PY'
import json, pprint
path='/Volumes/SSK SSD/TissueLab Revision/LFB/Seed48_Split/Training/autoresearch_runs/run_20260420_215745_90b105/round_0007/round_0007_worker/results.json'
obj=json.load(open(path))
for k in ['best_variation','feature_column','partial_r','selection_score','loo_predictive_r','is_loo_gap','adjusted_score','ranked_variations']:
    print(f'{k}:')
    print(json.dumps(obj.get(k), indent=2)[:8000])
PY
echo
echo '--- shared stats signatures ---'
python - <<'PY'
import inspect
from shared_analysis import stats, artifacts, sea_ad_lfb
for mod, names in [(stats,['partial_correlation','leave_one_out_summary','residualized_loo_predictive_correlation']),
                   (artifacts,['write_donor_feature_table','build_results_payload','write_results_payload']),
                   (sea_ad_lfb,['load_training_cohort','build_cell_table','load_region_polygons'])]:
    print(f'\nMODULE {mod.__name__}')
    for name in names:
        fn=getattr(mod,name)
        print(f'\n{name}{inspect.signature(fn)}')
        ds=inspect.getdoc(fn) or ''
        print(ds[:700])
PY
