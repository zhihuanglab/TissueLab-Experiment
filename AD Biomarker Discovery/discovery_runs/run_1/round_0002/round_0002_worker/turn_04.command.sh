set -euo pipefail
python - <<'PY'
import sys, inspect
sys.path.insert(0,'/shared/lib')
from shared_analysis import build_cell_table, load_region_polygons
print('build_cell_table signature:', inspect.signature(build_cell_table))
src = inspect.getsource(build_cell_table)
print(src[:5000])
print('--- load_region_polygons signature:', inspect.signature(load_region_polygons))
print(inspect.getsource(load_region_polygons)[:5000])
PY
echo '--- inspect one zarr tree shallow ---'
python - <<'PY'
from pathlib import Path
import zarr
slide=Path('/data/H19.33.004-A12-LFB.svs.zarr')
root=zarr.open_group(str(slide), mode='r')
def walk(g,prefix='',depth=0,maxdepth=3):
    if depth>maxdepth: return
    for k,v in g.groups():
        print(prefix+k+'/')
        walk(v,prefix+'  ',depth+1,maxdepth)
    for k,v in g.arrays():
        print(prefix+k, v.shape, v.dtype)
walk(root)
PY
