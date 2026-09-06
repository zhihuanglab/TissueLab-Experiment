set -euo pipefail
echo '--- result.py header ---'
sed -n '1,40p' /scratch/result.py
echo '--- results.json ---'
python - <<'PY'
import json, pathlib
for p in ['/scratch/results.json','/scratch/report.md','/scratch/donor_feature_table.csv']:
    path=pathlib.Path(p)
    print(p, path.exists(), path.stat().st_size if path.exists() else None)
if pathlib.Path('/scratch/results.json').exists():
    obj=json.loads(pathlib.Path('/scratch/results.json').read_text())
    print('best_variation', obj.get('best_variation'))
    print('feature_column', obj.get('feature_column'))
    print('ranked_n', len(obj.get('ranked_variations', [])))
PY
