set -euo pipefail
python - <<'PY'
import sys, os
sys.path.insert(0,'/shared/lib')
from pathlib import Path
from shared_analysis.sea_ad_lfb import build_cell_table, load_training_cohort
import numpy as np
data_root=Path('/data')
row=load_training_cohort(data_root).iloc[0]
cells=build_cell_table(data_root/row['slide_name'], include_regions=False)
print('x range', float(cells.x.min()), float(cells.x.max()))
print('y range', float(cells.y.min()), float(cells.y.max()))
try:
    import openslide
    slide=openslide.OpenSlide(str(data_root/row['slide_name'].replace('.zarr','')))
    print('slide dims', slide.dimensions)
    print('mpp', slide.properties.get('aperio.MPP'))
except Exception as e:
    print('openslide failed', repr(e))
PY
