set -euo pipefail
echo '--- result.py head ---'
sed -n '1,40p' /scratch/result.py
echo '--- results.json ---'
python - <<'PY'
import json, pathlib
p=pathlib.Path('/scratch/results.json')
obj=json.loads(p.read_text())
print(json.dumps(obj, indent=2)[:12000])
PY
echo '--- report.md ---'
sed -n '1,240p' /scratch/report.md
