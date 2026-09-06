set -euo pipefail
python - <<'PY'
from pathlib import Path
import inspect
from shared_analysis import build_cell_table
import shared_analysis
from shared_analysis import sea_ad_lfb
print('shared_analysis package:', shared_analysis.__file__)
print('build_cell_table signature:', inspect.signature(build_cell_table))
print('\nsource build_cell_table:\n')
src=inspect.getsource(build_cell_table)
print(src[:6000])
PY
echo '--- previous round 5 result.py ---'
sed -n '1,260p' /data/autoresearch_runs/run_20260420_215745_90b105/round_0005/round_0005_worker/result.py || true
echo '--- previous round 5 report.md ---'
sed -n '1,220p' /data/autoresearch_runs/run_20260420_215745_90b105/round_0005/round_0005_worker/report.md || true
