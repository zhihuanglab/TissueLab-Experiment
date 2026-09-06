set -euo pipefail
echo '--- context_bundle.json ---'
python - <<'PY'
import json, pprint
with open('/scratch/context_bundle.json') as f:
    obj=json.load(f)
pprint.pp(obj)
PY
echo '--- program.md ---'
sed -n '1,220p' /scratch/program.md
echo '--- existing result.py template ---'
sed -n '1,260p' /scratch/result.py || true
