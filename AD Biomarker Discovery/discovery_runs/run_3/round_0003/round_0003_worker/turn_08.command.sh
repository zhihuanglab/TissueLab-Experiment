set -euo pipefail
python - <<'PY'
import sys, inspect, json
from pathlib import Path
sys.path.append('/shared/lib')
from shared_analysis.stats import partial_correlation, residualized_loo_predictive_correlation
from shared_analysis.artifacts import write_donor_feature_table, write_results_payload, build_results_payload
print('partial_correlation sig', inspect.signature(partial_correlation))
print(inspect.getsource(partial_correlation)[:2500])
print('\nresidualized_loo_predictive_correlation sig', inspect.signature(residualized_loo_predictive_correlation))
print(inspect.getsource(residualized_loo_predictive_correlation)[:4500])
print('\nwrite_donor_feature_table sig', inspect.signature(write_donor_feature_table))
print(inspect.getsource(write_donor_feature_table)[:1800])
print('\nslide_manifest top-level type preview:')
manifest=json.load(open('/shared/cache/slide_manifest.json'))
print(type(manifest))
print(list(manifest.keys())[:10])
print('slide_count', manifest.get('slide_count'))
print('first slide keys', manifest['slides'][0].keys())
print('first slide mpp', manifest['slides'][0].get('mpp'))
PY
