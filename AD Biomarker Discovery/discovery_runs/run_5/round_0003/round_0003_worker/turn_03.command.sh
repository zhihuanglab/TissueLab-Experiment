set -euo pipefail
python - <<'PY'
import inspect, textwrap, json
from shared_analysis import build_cell_table, load_training_cohort
from shared_analysis import sea_ad_lfb, stats
print('build_cell_table signature:', inspect.signature(build_cell_table))
src = inspect.getsource(build_cell_table)
print('--- build_cell_table source excerpt ---')
print('\n'.join(src.splitlines()[:120]))
print('--- load_training_cohort signature:', inspect.signature(load_training_cohort))
print('--- stats funcs ---')
for name in ['partial_correlation','leave_one_out_summary','residualized_loo_predictive_correlation']:
    fn=getattr(stats,name)
    print(name, inspect.signature(fn))
    print('\n'.join(inspect.getsource(fn).splitlines()[:80]))
    print('-----')
PY
echo '--- quick sample cell table columns ---'
python - <<'PY'
from pathlib import Path
from shared_analysis import load_training_cohort, build_cell_table
cohort=load_training_cohort('/data')
slide=cohort.iloc[0]['slide_name']
cells=build_cell_table(Path('/data')/slide, include_regions=True, include_geometry=False)
print('slide', slide, 'shape', cells.shape)
print(cells.columns.tolist())
print(cells[['cell_label','region_label','centroid_x','centroid_y']].head().to_string())
print(cells['region_label'].value_counts(dropna=False).head(15).to_string())
print(cells['cell_label'].value_counts().to_string())
PY
