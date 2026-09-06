set -euo pipefail
echo '--- context_bundle.json ---'
python - <<'PY'
import json, textwrap
with open('/scratch/context_bundle.json') as f:
    obj=json.load(f)
import pprint
pprint.pp(obj)
PY
echo '--- program.md ---'
sed -n '1,240p' /scratch/program.md
echo '--- existing result.py template ---'
sed -n '1,240p' /scratch/result.py || true
