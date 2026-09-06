set -e
python - <<'PY'
from pathlib import Path
import pandas as pd, zarr, json, numpy as np
cohort = pd.read_csv('/data/training_cohort.csv')
slide = Path('/data') / cohort.iloc[0]['slide_name']
root = zarr.open_group(str(slide), mode='r')
print("centroids shape", root['SegmentationNode']['centroids'].shape, root['SegmentationNode']['centroids'].dtype)
print("class_id shape", root['ClassificationNode']['nuclei_class_id'].shape, root['ClassificationNode']['nuclei_class_id'].dtype)
print("class_name shape", root['ClassificationNode']['nuclei_class_name'].shape, root['ClassificationNode']['nuclei_class_name'].dtype)
print("first class names", root['ClassificationNode']['nuclei_class_name'][:10].tolist())
print("first class ids", root['ClassificationNode']['nuclei_class_id'][:10].tolist())
print("metadata keys", list(root['ClassificationNode']['metadata'].attrs.keys()))
for ann in ['CA1_01','CA2_01']:
    arr=root['CustomAnnotations'][ann]
    print(ann, "shape", arr.shape, "dtype", arr.dtype, "head", arr[:5].tolist())
PY
