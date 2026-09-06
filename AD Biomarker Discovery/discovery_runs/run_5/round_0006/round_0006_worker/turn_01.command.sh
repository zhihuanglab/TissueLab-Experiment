set -euo pipefail
echo '--- context_bundle.json ---'
python - <<'PY'
import json, pprint
with open('/scratch/context_bundle.json') as f:
    data=json.load(f)
pprint.pp(data)
PY
echo '--- program.md ---'
sed -n '1,240p' /scratch/program.md
echo '--- result.py template head ---'
sed -n '1,240p' /scratch/result.py || true
