set -euo pipefail
python - <<'PY'
import inspect, textwrap
from shared_analysis import build_cell_table, load_training_cohort
from shared_analysis import stats as st
import shared_analysis.sea_ad_lfb as sa
print('build_cell_table signature:', inspect.signature(build_cell_table))
print(textwrap.dedent(inspect.getsource(build_cell_table))[:5000])
print('\npartial_correlation sig:', inspect.signature(st.partial_correlation))
print(textwrap.dedent(inspect.getsource(st.partial_correlation))[:2000])
print('\nleave_one_out_summary sig:', inspect.signature(st.leave_one_out_summary))
print(textwrap.dedent(inspect.getsource(st.leave_one_out_summary))[:2600])
PY
python - <<'PY'
from shared_analysis import build_cell_table
from pathlib import Path
slide=Path('/data')/'H19.33.004-A12-LFB.svs.zarr'
df=build_cell_table(slide, include_regions=True, include_geometry=False)
print(df.head().to_string())
print(df.columns.tolist())
print(df[['cell_label','region_name']].value_counts().head(20).to_string())
PY
