set -euo pipefail
python - <<'PY'
from pathlib import Path
from shared_analysis import build_cell_table
slide=Path('/data/H19.33.004-A12-LFB.svs.zarr')
cells=build_cell_table(slide, include_regions=True, include_geometry=False)
print('columns:', list(cells.columns))
print(cells[['label','x_centroid','y_centroid','region_label']].head().to_string())
print('n', len(cells), 'labels', cells['label'].value_counts().to_dict())
print('regions sample', cells['region_label'].value_counts().head(10).to_dict())
PY
echo '--- round8 result.py ---'
sed -n '1,260p' "/data/autoresearch_runs/run_20260421_014751_94b13b/round_0008/round_0008_worker/result.py" || true
echo '--- round8 results.json ---'
python - <<'PY'
import json, pathlib
p=pathlib.Path('/data/autoresearch_runs/run_20260421_014751_94b13b/round_0008/round_0008_worker/results.json')
obj=json.loads(p.read_text())
print(json.dumps(obj, indent=2)[:12000])
PY
