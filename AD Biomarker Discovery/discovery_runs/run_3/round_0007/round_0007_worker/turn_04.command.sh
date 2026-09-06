set -euo pipefail
python - <<'PY'
from pathlib import Path
import inspect
import numpy as np
from shared_analysis import build_cell_table, load_centroids, load_contours, load_region_polygons, load_class_lookup
slide = Path('/data/H19.33.004-A12-LFB.svs.zarr')
print('--- build_cell_table sample ---')
cells = build_cell_table(slide, include_regions=True, include_geometry=False)
print(cells.columns.tolist())
print(cells[['x','y','cell_type','region']].head().to_string(index=False))
print('n=', len(cells))
print('\n--- low-level shapes ---')
print('centroids', load_centroids(slide).shape)
print('contours', load_contours(slide).shape)
polys = load_region_polygons(slide)
print('regions', sorted(polys.keys()))
print('CA1 polygon count', len(polys.get('CA1', [])))
lookup = load_class_lookup(slide)
print('class lookup', lookup)
print('\n--- leave_one_out_summary signature ---')
from shared_analysis import leave_one_out_summary, partial_correlation, residualized_loo_predictive_correlation, bootstrap_partial_correlation
for fn in [partial_correlation, leave_one_out_summary, residualized_loo_predictive_correlation, bootstrap_partial_correlation]:
    print(fn.__name__, inspect.signature(fn))
    doc = inspect.getdoc(fn)
    if doc:
        print(doc.splitlines()[0])
PY
