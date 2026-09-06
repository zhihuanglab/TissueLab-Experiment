set -euo pipefail
python - <<'PY'
import json, glob, os, math
for f in glob.glob('/data/autoresearch_runs/**/results.json', recursive=True)[:200]:
    try:
        d=json.load(open(f))
    except: 
        continue
    rv=d.get('ranked_variations')
    if not rv: continue
    for row in rv[:3]:
        if all(k in row for k in ['partial_r','loo_predictive_r','gap','gap_penalty','adjusted_score']):
            pr=abs(row['partial_r']); loo=abs(row['loo_predictive_r']); gap=row['gap']; pen=row['gap_penalty']; adj=row['adjusted_score']
            print(os.path.basename(os.path.dirname(f)), row['name'], 'pr',pr,'loo',loo,'gap',gap,'pen',pen,'adj',adj, 'adj-loo', adj-loo)
            break
    break
PY
