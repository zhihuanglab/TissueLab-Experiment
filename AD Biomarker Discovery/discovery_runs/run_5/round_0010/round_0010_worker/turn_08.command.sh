set -euo pipefail
python - <<'PY'
import sys, inspect, json, glob
sys.path.insert(0,'/shared/lib')
from shared_analysis import stats
for name in ['partial_correlation','leave_one_out_summary','residualized_loo_predictive_correlation']:
    fn=getattr(stats,name)
    print(f'\n{name}{inspect.signature(fn)}')
    src=inspect.getsource(fn)
    print('\n'.join(src.splitlines()[:80]))
print('\n--- sample previous round results keys ---')
import pathlib
for p in glob.glob('/Volumes/SSK SSD/TissueLab Revision/LFB/Seed48_Split/Training/autoresearch_runs/run_20260421_052101_efbee1/round_0009/round_0009_worker/results.json'):
    data=json.load(open(p))
    print(pathlib.Path(p).name, data.keys())
    print(json.dumps({k:data[k] for k in data if k not in ['per_donor_loo','ranked_variations']}, indent=2)[:2000])
PY
