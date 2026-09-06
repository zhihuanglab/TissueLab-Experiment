set -euo pipefail
python - <<'PY'
import sys, json, pathlib
sys.path.insert(0,'/shared/lib')
from pathlib import Path
from shared_analysis.sea_ad_lfb import build_cell_table, load_training_cohort
from shared_analysis.stats import partial_correlation, residualized_loo_predictive_correlation
import inspect

print('--- stats helper signatures ---')
import shared_analysis.stats as st
for name in ['partial_correlation','residualized_loo_predictive_correlation','leave_one_out_summary']:
    fn=getattr(st,name)
    print(name, inspect.signature(fn))
    src=inspect.getsource(fn).splitlines()
    print('\n'.join(src[:40]))
    print('...')

print('\n--- artifacts helpers ---')
import shared_analysis.artifacts as art
for name in ['write_donor_feature_table','build_results_payload','write_results_payload','coerce_results_payload']:
    fn=getattr(art,name)
    print(name, inspect.signature(fn))

print('\n--- nearby prior result files under /data/autoresearch_runs ---')
base=Path('/data/autoresearch_runs')
matches=list(base.rglob('round_0006_worker/result.py'))
print('matches:', len(matches))
for p in matches[:3]:
    print(p)
    print('--- head ---')
    print('\n'.join(p.read_text().splitlines()[:220]))
PY
