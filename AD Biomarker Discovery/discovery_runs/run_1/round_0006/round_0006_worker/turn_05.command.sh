set -euo pipefail
echo '--- previous round scripts/results snippets ---'
python - <<'PY'
from pathlib import Path
import json, glob
paths=sorted(glob.glob('/data/autoresearch_runs/*/round_0004/round_0004_worker/results.json'))
print('found', len(paths))
if paths:
    p=Path(paths[-1]); print('using', p)
    print(p.read_text()[:4000])
    rp=p.with_name('result.py')
    print('\n--- result.py head ---')
    print(rp.read_text()[:5000])
PY

echo '--- stats/artifacts snippets ---'
python - <<'PY'
from pathlib import Path
for path in ['/shared/lib/shared_analysis/stats.py','/shared/lib/shared_analysis/artifacts.py']:
    text=Path(path).read_text()
    print(f'\n### {path}')
    for needle in ['def partial_correlation','def residualized_loo_predictive_correlation','def write_donor_feature_table','def build_results_payload','def validate_results_payload']:
        idx=text.find(needle)
        if idx!=-1:
            print(text[idx:idx+2200])
PY
