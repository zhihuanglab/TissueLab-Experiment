set -euo pipefail
python - <<'PY'
import sys
sys.path.append('/shared/lib')
import zarr, pandas as pd
cohort = pd.read_csv('/data/training_cohort.csv')
path = '/data/' + cohort.iloc[0]['slide_name']
root = zarr.open(path, mode='r')
print(root['SegmentationNode']['centroids'].shape, root['SegmentationNode']['centroids'].dtype)
print(root['SegmentationNode']['contours'].shape, root['SegmentationNode']['contours'].dtype)
print(root['ClassificationNode']['nuclei_class_id'].shape, root['ClassificationNode']['nuclei_class_id'].dtype)
PY
