set -euo pipefail
python - <<'PY'
import json, pandas as pd, os, sys
from pathlib import Path
sys.path.append('/shared/lib')
from shared_analysis.sea_ad_lfb import load_class_lookup, load_centroids, load_class_ids, load_region_annotations

cohort = pd.read_csv('/data/training_cohort.csv')
row = cohort.iloc[0]
slide = Path('/data') / row['slide_name']
print('sample slide:', slide)
print('columns:', cohort.columns.tolist())
print('shape', cohort.shape)

print('\nslide manifest sample:')
manifest = json.load(open('/shared/cache/slide_manifest.json'))
for k,v in list(manifest.items())[:1]:
    print(k, v)
    break

print('\nclass lookup:', load_class_lookup(slide))
cent = load_centroids(slide)
cid = load_class_ids(slide)
print('centroids shape', cent.shape, cent.dtype, 'class_ids', cid.shape, cid.dtype)
anns = load_region_annotations(slide, scale=16.0)
print('first 5 annotations:')
for a in anns[:5]:
    print({kk:(vv[:2].tolist() if kk=="points" else vv) for kk,vv in a.items()})

print('\nshared cached region polygon file snippet:')
cache_path = Path('/shared/cache/region_polygons') / f"{row['slide_name']}.json"
print(cache_path)
print(open(cache_path).read()[:500])
PY
