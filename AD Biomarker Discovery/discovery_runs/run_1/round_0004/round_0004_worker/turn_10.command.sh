set -euo pipefail
echo '--- build_cell_table source excerpt ---'
sed -n '1,260p' /shared/lib/shared_analysis/sea_ad_lfb.py | sed -n '1,240p'
echo
echo '--- stats/artifacts source excerpts ---'
grep -n "def partial_correlation\|def leave_one_out_summary\|def residualized_loo_predictive_correlation" -n /shared/lib/shared_analysis/stats.py
grep -n "def write_donor_feature_table\|def build_results_payload\|def write_results_payload\|def validate_results_payload" -n /shared/lib/shared_analysis/artifacts.py
echo
cat > /scratch/_sample.py <<'PY'
from pathlib import Path
from shared_analysis import load_training_cohort, build_cell_table
cohort = load_training_cohort('/data')
slide = cohort.iloc[0]['slide_name']
cells = build_cell_table(Path('/data')/slide, include_regions=True, include_geometry=False)
print('slide', slide, 'shape', cells.shape)
print('columns', cells.columns.tolist())
print(cells[['cell_type','region','x_centroid','y_centroid']].head().to_string())
print('region counts top:')
print(cells['region'].value_counts(dropna=False).head(12).to_string())
print('cell_type counts top:')
print(cells['cell_type'].value_counts().head(12).to_string())
print('CA1 types top:')
ca1 = cells[cells['region']=='CA1']
print(ca1['cell_type'].value_counts().head(12).to_string())
PY
PYTHONPATH=/shared/lib /usr/local/bin/python3 /scratch/_sample.py
