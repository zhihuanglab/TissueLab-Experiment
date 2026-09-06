set -e
python - <<'PY'
import json, glob, os
for p in sorted(glob.glob('/data/autoresearch_runs/*/round_0004/round_0004_worker/results.json'))[:3]:
    print("FILE", p)
    obj=json.load(open(p))
    print(json.dumps({k: obj.get(k) for k in ['best_variation','partial_r','selection_score','loo_predictive_r','is_loo_gap','adjusted_score','feature_column']}, indent=2))
    print("rank head", obj.get('ranked_variations',[])[:3])
PY
