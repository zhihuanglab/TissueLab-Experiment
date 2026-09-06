set -euo pipefail
python - <<'PY'
import sys, inspect, textwrap, json
sys.path.append('/shared/lib')
import shared_analysis.sea_ad_lfb as sa
import shared_analysis.stats as st
for name in ['load_training_cohort','build_slide_manifest','load_region_polygons','assign_centroids_to_regions','load_centroids','load_class_ids','load_class_lookup']:
    fn=getattr(sa,name)
    print(f'\n## {name} {inspect.signature(fn)}')
    src=textwrap.dedent(inspect.getsource(fn))
    print(src[:1800])
for name in ['residualized_loo_predictive_correlation','leave_one_out_partial_correlation','residualize']:
    fn=getattr(st,name)
    print(f'\n## {name} {inspect.signature(fn)}')
    src=textwrap.dedent(inspect.getsource(fn))
    print(src[:2400])
PY
python - <<'PY'
import sys
sys.path.append('/shared/lib')
from shared_analysis import build_cell_table
from pathlib import Path
slide=Path('/data')/'H19.33.004-A12-LFB.svs.zarr'
df=build_cell_table(slide, include_regions=True, include_geometry=False)
print('rows', len(df), 'region non-null', df['region'].notna().sum())
print(df[['cell_type','region']].dropna().value_counts().head(20).to_string())
print('unique regions:', sorted(df['region'].dropna().unique().tolist()))
PY
