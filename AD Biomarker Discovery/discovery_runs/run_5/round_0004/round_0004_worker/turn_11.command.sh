set -euo pipefail
python - <<'PY'
import openslide, pandas as pd
cohort = pd.read_csv('/data/training_cohort.csv')
svs='/data/'+cohort.iloc[0]['slide_name'].replace('.zarr','')
slide=openslide.OpenSlide(svs)
print(slide.properties.get('openslide.mpp-x'), slide.properties.get('openslide.mpp-y'))
print([k for k in slide.properties.keys() if 'mpp' in k.lower()][:10])
PY
