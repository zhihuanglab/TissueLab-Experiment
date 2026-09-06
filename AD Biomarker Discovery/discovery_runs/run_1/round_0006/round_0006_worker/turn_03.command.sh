set -euo pipefail
echo '--- sample zarr tree (top 3 levels) ---'
python - <<'PY'
from pathlib import Path
import zarr
sample = Path('/data') / 'H20.33.001-A12-LFB.svs.zarr'
root = zarr.open(str(sample), mode='r')
def walk(group, prefix='', depth=0, maxdepth=3):
    if depth>maxdepth: return
    for k, v in group.groups():
        print(prefix + k + '/')
        walk(v, prefix + '  ', depth+1, maxdepth)
    for k, v in group.arrays():
        print(prefix + k, v.shape, v.dtype)
walk(root)
PY

echo '--- shared lib relevant defs ---'
python - <<'PY'
from pathlib import Path
for path in [
    Path('/shared/lib/shared_analysis/sea_ad_lfb.py'),
    Path('/shared/lib/shared_analysis/stats.py'),
    Path('/shared/lib/shared_analysis/artifacts.py'),
]:
    print(f'\n### {path}')
    text=path.read_text()
    for needle in ['def load_training_cohort','def build_cell_table','def load_region_polygons','def compute_contour_geometry',
                   'def partial_correlation','def leave_one_out_summary','def residualized_loo_predictive_correlation',
                   'def write_donor_feature_table']:
        idx=text.find(needle)
        if idx!=-1:
            print(text[idx:idx+1800])
PY
