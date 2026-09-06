set -euo pipefail
python - <<'PY'
import json, glob, os
count=0
for f in glob.glob('/data/autoresearch_runs/**/results.json', recursive=True):
    try:
        d=json.load(open(f))
    except: 
        continue
    rv=d.get('ranked_variations')
    if not rv: continue
    for row in rv:
        keys=['partial_r','loo_predictive_r','gap','gap_penalty','adjusted_score']
        if all(k in row and row[k] is not None for k in keys):
            pr=abs(row['partial_r']); loo=abs(row['loo_predictive_r']); gap=row['gap']; pen=row['gap_penalty']; adj=row['adjusted_score']
            if pen and pen>0:
                print(f,row['name'],'pr',pr,'loo',loo,'gap',gap,'pen',pen,'adj',adj,'adj-vs-loo',adj-loo,'gap-0.15',gap-0.15,'gap*0.05',gap*0.05)
                count+=1
                if count>=8: raise SystemExit
PY
