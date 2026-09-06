PYTHONPATH=/shared/lib python - <<'PY'
import inspect, textwrap
from pathlib import Path
import zarr
from shared_analysis import build_cell_table, load_region_polygons
from shared_analysis.sea_ad_lfb import compute_contour_geometry

slide = Path('/data') / 'H19.33.004-A12-LFB.svs.zarr'
root = zarr.open(str(slide), mode='r')
print('top groups:', list(root.group_keys()))
for g in ['SegmentationNode','ClassificationNode','CustomAnnotations']:
    grp=root[g]
    print('\nGROUP', g)
    print('  subgroups:', list(grp.group_keys()))
    print('  arrays:', list(grp.array_keys())[:20])

print('\nload_region_polygons signature:', inspect.signature(load_region_polygons))
print('build_cell_table signature:', inspect.signature(build_cell_table))
print('compute_contour_geometry signature:', inspect.signature(compute_contour_geometry))
print('\nsource load_region_polygons:\n')
print(textwrap.dedent(inspect.getsource(load_region_polygons))[:2200])
print('\nsource build_cell_table head:\n')
print(textwrap.dedent(inspect.getsource(build_cell_table))[:3200])
PY
