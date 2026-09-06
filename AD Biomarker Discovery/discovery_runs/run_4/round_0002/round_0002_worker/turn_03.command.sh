set -e
python - <<'PY'
from pathlib import Path
import zarr, json, inspect
from shared_analysis import sea_ad_lfb
slide = Path('/data/H19.33.004-A12-LFB.svs.zarr')
print('--- zarr tree top 3 levels ---')
root = zarr.open(str(slide), mode='r')
def walk(g, prefix='', depth=0, maxdepth=3):
    if depth>maxdepth: return
    for k,v in g.groups():
        print(prefix + k + '/')
        walk(v, prefix+'  ', depth+1, maxdepth)
    for k,v in g.arrays():
        try: shp=v.shape; dt=v.dtype
        except Exception: shp='?'; dt='?'
        print(prefix + f'{k}  shape={shp} dtype={dt}')
walk(root)

print('\n--- shared_analysis.sea_ad_lfb helper signatures/source snippets ---')
for name in ['build_cell_table','load_region_polygons']:
    fn = getattr(sea_ad_lfb, name)
    print(f'\n{name}{inspect.signature(fn)}')
    src = inspect.getsource(fn)
    print(src[:3000])
PY
