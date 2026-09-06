set -euo pipefail
export PYTHONPATH=/shared/lib:${PYTHONPATH:-}
python - <<'PY'
import inspect
from shared_analysis import build_cell_table, load_region_polygons, compute_contour_geometry
print('=== build_cell_table ===')
print(inspect.getsource(build_cell_table))
print('\n=== load_region_polygons ===')
print(inspect.getsource(load_region_polygons))
print('\n=== compute_contour_geometry ===')
print(inspect.getsource(compute_contour_geometry))
PY
