set -euo pipefail
python - <<'PY'
import inspect, textwrap, json, os
from pathlib import Path
import zarr
from shared_analysis import build_cell_table, load_region_polygons
import shared_analysis.sea_ad_lfb as mod

slide = Path('/data') / 'H19.33.004-A12-LFB.svs.zarr'
print('exists', slide.exists())
root = zarr.open(str(slide), mode='r')
print('top groups:', list(root.group_keys()))
for g in ['SegmentationNode','ClassificationNode','CustomAnnotations']:
    if g in root:
        grp=root[g]
        print('\nGROUP', g, 'subkeys', list(grp.group_keys()), 'arrays', list(grp.array_keys())[:20])

print('\nload_region_polygons signature:', inspect.signature(load_region_polygons))
print('build_cell_table signature:', inspect.signature(build_cell_table))
print('\nsource load_region_polygons:')
print(textwrap.dedent(inspect.getsource(load_region_polygons))[:1800])
print('\nsource build_cell_table head:')
print(textwrap.dedent(inspect.getsource(build_cell_table))[:2600])
PY
