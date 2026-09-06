set -euo pipefail
python - <<'PY'
import inspect
from shared_analysis import build_cell_table, load_region_polygons, load_training_cohort, build_slide_manifest
print('build_cell_table signature:', inspect.signature(build_cell_table))
print(inspect.getsource(build_cell_table)[:4000])
print('\nload_region_polygons signature:', inspect.signature(load_region_polygons))
print(inspect.getsource(load_region_polygons)[:2500])
print('\nload_training_cohort signature:', inspect.signature(load_training_cohort))
print(inspect.getsource(load_training_cohort)[:1500])
print('\nbuild_slide_manifest signature:', inspect.signature(build_slide_manifest))
print(inspect.getsource(build_slide_manifest)[:2000])
PY
