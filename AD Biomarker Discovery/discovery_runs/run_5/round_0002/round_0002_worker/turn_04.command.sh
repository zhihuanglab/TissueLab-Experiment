set -euo pipefail
python - <<'PY'
import sys, inspect, textwrap
sys.path.insert(0, '/shared/lib')
from shared_analysis import build_cell_table, load_region_polygons
print('build_cell_table source:')
print(textwrap.dedent(inspect.getsource(build_cell_table)))
print('\nload_region_polygons source:')
print(textwrap.dedent(inspect.getsource(load_region_polygons)))
PY

echo '--- sample zarr tree and arrays ---'
python - <<'PY'
from pathlib import Path
import zarr
slide=Path('/data/H19.33.004-A12-LFB.svs.zarr')
root=zarr.open(str(slide), mode='r')
def walk(g,p=''):
    if hasattr(g,'group_keys'):
        for k in g.group_keys():
            print('GROUP', p+'/'+k)
            walk(g[k], p+'/'+k)
        for k in g.array_keys():
            arr=g[k]
            print('ARRAY', p+'/'+k, 'shape=', arr.shape, 'dtype=', arr.dtype)
walk(root)
PY
