set -e
python - <<'PY'
from pathlib import Path
for p in ['/scratch/result.py','/scratch/report.md','/scratch/results.json','/scratch/donor_feature_table.csv']:
    print(p, Path(p).exists(), Path(p).stat().st_size if Path(p).exists() else None)
PY
