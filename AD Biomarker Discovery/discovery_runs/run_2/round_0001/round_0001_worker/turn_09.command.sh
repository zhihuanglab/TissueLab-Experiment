python - <<'PY'
import sys, inspect, textwrap
sys.path.insert(0, '/shared/lib')
from shared_analysis.sea_ad_lfb import (
    load_region_annotations, assign_centroids_to_regions, load_training_cohort,
    load_class_lookup, load_class_ids, load_centroids
)
from shared_analysis import partial_correlation, leave_one_out_summary, residualized_loo_predictive_correlation
import shared_analysis.stats as stats
import shared_analysis.artifacts as artifacts

for fn in [load_region_annotations, assign_centroids_to_regions, load_class_lookup, load_training_cohort]:
    print(f'\n--- {fn.__name__} ---')
    print(textwrap.dedent(inspect.getsource(fn))[:2600])

print('\n--- stats.partial_correlation ---')
print(textwrap.dedent(inspect.getsource(partial_correlation))[:2200])
print('\n--- stats.residualized_loo_predictive_correlation ---')
print(textwrap.dedent(inspect.getsource(residualized_loo_predictive_correlation))[:2600])
print('\n--- stats.leave_one_out_summary ---')
print(textwrap.dedent(inspect.getsource(leave_one_out_summary))[:2600])
print('\n--- artifacts write donor and results payload fns ---')
for name in ['write_donor_feature_table','build_results_payload','write_results_payload']:
    fn = getattr(artifacts,name)
    print(f'\n{name}:')
    print(textwrap.dedent(inspect.getsource(fn))[:2200])
PY
