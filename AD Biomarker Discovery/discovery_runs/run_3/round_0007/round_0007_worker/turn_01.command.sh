set -euo pipefail
echo '--- context_bundle.json ---'
python - <<'PY'
import json
with open('/scratch/context_bundle.json') as f:
    obj=json.load(f)
print(json.dumps(obj, indent=2)[:20000])
PY
echo
echo '--- program.md ---'
sed -n '1,220p' /scratch/program.md
echo
echo '--- existing result.py template ---'
sed -n '1,260p' /scratch/result.py
