set -euo pipefail
echo '--- helper signatures / relevant implementations ---'
python - <<'PY'
import inspect
from shared_analysis import build_cell_table, load_region_polygons, compute_contour_geometry
from shared_analysis import partial_correlation, residualized_loo_predictive_correlation, leave_one_out_summary
for fn in [build_cell_table, load_region_polygons, compute_contour_geometry, partial_correlation, residualized_loo_predictive_correlation, leave_one_out_summary]:
    print(f"\n## {fn.__name__}")
    try:
        print(inspect.signature(fn))
        src = inspect.getsource(fn)
        print("\n".join(src.splitlines()[:80]))
    except Exception as e:
        print("could not inspect", e)
PY

echo '--- one-slide data sanity ---'
python - <<'PY'
from pathlib import Path
import pandas as pd
from shared_analysis import build_cell_table, load_region_polygons
cohort = pd.read_csv('/data/training_cohort.csv')
slide = Path('/data') / cohort.iloc[0]['slide_name']
print('slide', slide)
polys = load_region_polygons(slide)
print('region keys', list(polys.keys())[:10], 'n=', len(polys))
cells = build_cell_table(slide, include_regions=True, include_geometry=True)
print(cells.columns.tolist())
print(cells[['cell_label','region_label','centroid_x','centroid_y']].head().to_string())
print('regions value counts head')
print(cells['region_label'].value_counts(dropna=False).head(12).to_string())
print('geom cols', [c for c in cells.columns if 'contour' in c or 'area' in c or 'perimeter' in c][:20])
print(cells[[c for c in cells.columns if 'area' in c][:5]].head().to_string())
PY
