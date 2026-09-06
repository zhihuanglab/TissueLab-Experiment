set -euo pipefail
python - <<'PY'
import json, glob, math
shown=0
for p in glob.glob('/data/autoresearch_runs/**/results.json', recursive=True):
    try: obj=json.load(open(p))
    except: continue
    if 'selection_score' in obj and 'adjusted_score' in obj and 'is_loo_gap' in obj:
        vals={k:obj.get(k) for k in ['feature_name','selection_score','partial_r','loo_predictive_r','is_loo_gap','adjusted_score']}
        if all(isinstance(v,(int,float)) for k,v in vals.items() if k!='feature_name'):
            print(p)
            print(vals)
            shown+=1
            if shown>=5: break
PY
