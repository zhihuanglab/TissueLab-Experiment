set -euo pipefail
python - <<'PY'
import sys, inspect
sys.path.append('/shared/lib')
from shared_analysis import build_cell_table, load_region_polygons, load_training_cohort
print('build_cell_table signature:', inspect.signature(build_cell_table))
print(inspect.getsource(build_cell_table)[:5000])
print('\nload_region_polygons signature:', inspect.signature(load_region_polygons))
print(inspect.getsource(load_region_polygons)[:3000])
PY
