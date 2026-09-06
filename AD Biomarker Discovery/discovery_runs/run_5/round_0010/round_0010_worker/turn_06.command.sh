set -euo pipefail
python - <<'PY'
import sys
sys.path.insert(0, '/shared/lib')
import inspect, zarr
from pathlib import Path
from shared_analysis import sea_ad_lfb
for name in ['build_cell_table','load_training_cohort','load_region_polygons','compute_contour_geometry']:
    fn=getattr(sea_ad_lfb,name)
    print(f'\n{name}{inspect.signature(fn)}')
    src=inspect.getsource(fn)
    print('\n'.join(src.splitlines()[:55]))
sample=Path('/data/H19.33.004-A12-LFB.svs.zarr')
print('\n--- group keys top ---')
root=zarr.open(str(sample), mode='r')
print(list(root.group_keys()))
for g in root.group_keys():
    grp=root[g]
    print(g, 'arrays', list(grp.array_keys())[:10], 'groups', list(grp.group_keys())[:10])

print('\n--- build_cell_table sample columns ---')
cells=sea_ad_lfb.build_cell_table(sample, include_regions=True, include_geometry=True)
print(cells.columns.tolist())
print(cells[['x_centroid','y_centroid','cell_type','region_label']].head())
print(cells[['contour_area_pixels']].describe())
PY
echo '--- previous round 7/9 script snippets ---'
python - <<'PY'
from pathlib import Path
for p in [
Path('/Volumes/SSK SSD/TissueLab Revision/LFB/Seed48_Split/Training/autoresearch_runs/run_20260421_052101_efbee1/round_0007/round_0007_worker/result.py'),
Path('/Volumes/SSK SSD/TissueLab Revision/LFB/Seed48_Split/Training/autoresearch_runs/run_20260421_052101_efbee1/round_0009/round_0009_worker/result.py'),
]:
    txt=p.read_text()
    print(f"\n--- {p.name} {p.parent.name} ---")
    for key in ['def compute_donor_score', 'from scipy.spatial import cKDTree', 'FEATURE_NAME']:
        i=txt.find(key)
        if i>=0:
            print(txt[i:i+1800])
            break
PY
