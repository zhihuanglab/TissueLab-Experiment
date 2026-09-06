set -euo pipefail
python - <<'PY'
import sys, inspect
sys.path.append('/shared/lib')
import shared_analysis
print('shared_analysis', shared_analysis.__file__)
from shared_analysis import build_cell_table, load_training_cohort, stats
print('build_cell_table signature:', inspect.signature(build_cell_table))
print('--- build_cell_table excerpt ---')
print('\n'.join(inspect.getsource(build_cell_table).splitlines()[:140]))
cohort=load_training_cohort('/data')
slide=cohort.iloc[0]['slide_name']
cells=build_cell_table('/data/'+slide, include_regions=True, include_geometry=False)
print('slide', slide, 'shape', cells.shape)
print('cell columns', cells.columns.tolist())
print(cells[['cell_label','region_label','centroid_x','centroid_y']].head().to_string())
print('regions', cells['region_label'].value_counts(dropna=False).head(10).to_dict())
print('labels', cells['cell_label'].value_counts().to_dict())
print('partial_correlation sig', inspect.signature(stats.partial_correlation))
print('\n'.join(inspect.getsource(stats.partial_correlation).splitlines()[:120]))
print('loo sig', inspect.signature(stats.residualized_loo_predictive_correlation))
print('\n'.join(inspect.getsource(stats.residualized_loo_predictive_correlation).splitlines()[:150]))
PY
