set -e
python - <<'PY'
import inspect
from shared_analysis import build_cell_table, load_region_polygons, load_training_cohort
import shared_analysis
import shared_analysis.sea_ad_lfb as m
print('module path', m.__file__)
print('\nload_training_cohort:\n', inspect.getsource(load_training_cohort))
print('\nload_region_polygons signature', inspect.signature(load_region_polygons))
print('\nload_region_polygons head:\n', '\n'.join(inspect.getsource(load_region_polygons).splitlines()[:80]))
print('\nbuild_cell_table signature', inspect.signature(build_cell_table))
src=inspect.getsource(build_cell_table).splitlines()
print('\nbuild_cell_table head:')
print('\n'.join(src[:140]))
PY
