set -euo pipefail
python - <<'PY'
import sys, warnings
sys.path.insert(0, '/shared/lib')
from pathlib import Path
import numpy as np
from matplotlib.path import Path as MplPath
from shared_analysis import load_training_cohort
import shared_analysis.sea_ad_lfb as sea
warnings.filterwarnings('ignore')
cohort=load_training_cohort('/data')
slide_name=cohort.iloc[0]['slide_name']
z=Path('/data')/slide_name
cent=sea.load_centroids(z)
class_ids=sea.load_class_ids(z)
lookup=sea.load_class_lookup(z)
polys=sea.load_region_polygons(z)
ca1=polys.get('CA1', [])
mask=np.zeros(len(cent), dtype=bool)
for poly in ca1:
    mask |= MplPath(poly, closed=True).contains_points(cent)
unique, counts=np.unique(class_ids[mask], return_counts=True)
print('slide', slide_name, 'n_ca1', int(mask.sum()), 'n_polys', len(ca1))
print({lookup[int(u)]: int(c) for u,c in zip(unique,counts)})
PY
