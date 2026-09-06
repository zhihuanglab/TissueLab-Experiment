set -euo pipefail
python - <<'PY'
import sys, inspect
sys.path.append('/shared/lib')
from pathlib import Path
from shared_analysis.sea_ad_lfb import build_cell_table, load_training_cohort
from shared_analysis.stats import partial_correlation, leave_one_out_summary, residualized_loo_predictive_correlation
from shared_analysis.artifacts import write_donor_feature_table, write_results_payload

data_root = Path('/data')
cohort = load_training_cohort(data_root)
slide_name = cohort.iloc[0]['slide_name']
cells = build_cell_table(data_root/slide_name, include_regions=True, include_geometry=False)
print('cell_type unique sample:', sorted(cells['cell_type'].dropna().unique().tolist()))
print('region unique:', sorted(cells['region'].dropna().unique().tolist()))
print('CA1 pyramidal n:', int(((cells['region']=='CA1') & (cells['cell_type']=='Pyramidal Neuron')).sum()))
print('CA1 reactive astrocyte n:', int(((cells['region']=='CA1') & (cells['cell_type']=='Reactive Astrocyte')).sum()))
print('partial_correlation sig:', inspect.signature(partial_correlation))
print('leave_one_out_summary sig:', inspect.signature(leave_one_out_summary))
print('residualized_loo_predictive_correlation sig:', inspect.signature(residualized_loo_predictive_correlation))
print('write_donor_feature_table sig:', inspect.signature(write_donor_feature_table))
PY
