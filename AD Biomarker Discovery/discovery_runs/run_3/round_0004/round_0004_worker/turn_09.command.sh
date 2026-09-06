set -euo pipefail
python - <<'PY'
from pathlib import Path
import json, math
files = list(Path('/data/autoresearch_runs').rglob('results.json'))[:200]
seen=0
for p in files:
    try:
        d=json.loads(p.read_text())
    except Exception:
        continue
    if all(k in d for k in ['partial_r','loo_predictive_r','adjusted_score','is_loo_gap']):
        pr=d['partial_r']; loo=d['loo_predictive_r']; adj=d['adjusted_score']; gap=d['is_loo_gap']; pen=d.get('gap_penalty')
        if pr is not None and loo is not None and adj is not None and gap is not None:
            print(p)
            print('partial',pr,'loo',loo,'adj',adj,'gap',gap,'pen',pen)
            seen+=1
            if seen>=8: break
PY
