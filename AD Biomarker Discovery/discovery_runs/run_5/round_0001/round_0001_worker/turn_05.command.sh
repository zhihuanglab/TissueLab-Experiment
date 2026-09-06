set -euo pipefail
export PYTHONPATH=/shared/lib:${PYTHONPATH:-}
echo '--- shared_analysis.__init__ ---'
sed -n '1,220p' /shared/lib/shared_analysis/__init__.py
echo '--- build_cell_table def ---'
python - <<'PY'
import inspect
from shared_analysis.sea_ad_lfb import build_cell_table, load_training_cohort
print(inspect.signature(load_training_cohort))
print(inspect.signature(build_cell_table))
print(inspect.getsource(build_cell_table)[:5000])
PY
echo '--- sample cell table summary ---'
python - <<'PY'
from pathlib import Path
from shared_analysis.sea_ad_lfb import load_training_cohort, build_cell_table
data_root=Path('/data')
cohort=load_training_cohort(data_root)
slide = data_root / cohort.iloc[0]['slide_name']
cells = build_cell_table(slide, include_regions=True, include_geometry=False)
print(cells[['cell_type','region']].head(8).to_string())
print('cell types:', cells['cell_type'].dropna().astype(str).value_counts().head(12).to_string())
print('regions:', cells['region'].dropna().astype(str).value_counts().to_string())
PY
