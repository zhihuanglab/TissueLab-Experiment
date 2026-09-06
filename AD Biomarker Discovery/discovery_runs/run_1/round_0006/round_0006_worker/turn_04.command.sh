set -euo pipefail
echo '--- /data head ---'
find /data -maxdepth 1 \( -name '*.zarr' -o -name 'training_cohort.csv' \) | sed 's#/data/##' | sort | head -20

echo '--- first slide_name from cohort ---'
python - <<'PY'
import pandas as pd
df=pd.read_csv('/data/training_cohort.csv')
print(df[['donor_id','slide_name']].head(5).to_string(index=False))
PY

slide=$(python - <<'PY'
import pandas as pd
df=pd.read_csv('/data/training_cohort.csv')
print(df['slide_name'].iloc[0])
PY
)
echo "--- sample slide: $slide ---"
python - <<'PY'
from pathlib import Path
import pandas as pd, zarr
slide = pd.read_csv('/data/training_cohort.csv')['slide_name'].iloc[0]
root = zarr.open(str(Path('/data')/slide), mode='r')
def walk(group, prefix='', depth=0, maxdepth=2):
    if depth>maxdepth: return
    for k, v in group.groups():
        print(prefix + k + '/')
        walk(v, prefix + '  ', depth+1, maxdepth)
    for k, v in group.arrays():
        print(prefix + k, v.shape, v.dtype)
walk(root)
PY

echo '--- shared lib snippets ---'
python - <<'PY'
from pathlib import Path
path=Path('/shared/lib/shared_analysis/sea_ad_lfb.py')
text=path.read_text()
for needle in ['def load_training_cohort','def build_cell_table','def load_region_polygons','def compute_contour_geometry']:
    idx=text.find(needle)
    print(f'\n## {needle}')
    print(text[idx:idx+2200])
PY
