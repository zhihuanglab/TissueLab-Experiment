set -euo pipefail
python - <<'PY'
import sys, numpy as np, zarr
sys.path.insert(0,'/shared/lib')
from shared_analysis import load_training_cohort
cohort=load_training_cohort('/data')
slide='/data/'+cohort.iloc[0]['slide_name']
root=zarr.open(slide, mode='r')
for arr_name in ['SegmentationNode/centroids','SegmentationNode/contours','ClassificationNode/nuclei_class_id']:
    arr=root[arr_name]
    print(arr_name, 'shape=', arr.shape, 'dtype=', arr.dtype, 'chunks=', arr.chunks)
PY
