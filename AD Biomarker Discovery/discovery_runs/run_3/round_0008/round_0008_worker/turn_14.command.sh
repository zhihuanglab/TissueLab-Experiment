set -euo pipefail
python - <<'PY'
import json, glob, math
count=0
for p in glob.glob('/data/autoresearch_runs/**/results.json', recursive=True):
    try: obj=json.load(open(p))
    except: continue
    if 'selection_score' in obj and 'adjusted_score' in obj and 'is_loo_gap' in obj:
        ss=obj['selection_score']; gap=obj['is_loo_gap']; adj=obj['adjusted_score']
        if all(isinstance(x,(int,float)) for x in [ss,gap,adj]) and math.isfinite(ss) and math.isfinite(gap) and math.isfinite(adj):
            print(ss, gap, adj, 'ss-gap', ss-gap, 'ss-abs', ss-abs(gap))
            count+=1
            if count>=5: break
PY
