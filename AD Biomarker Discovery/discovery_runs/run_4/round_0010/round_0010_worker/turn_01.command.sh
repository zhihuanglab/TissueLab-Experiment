set -euo pipefail
echo '--- context_bundle.json ---'
python - <<'PY'
import json, pprint
with open('/scratch/context_bundle.json') as f:
    data=json.load(f)
pprint.pp(data)
PY
echo '--- program.md ---'
cat /scratch/program.md
echo '--- existing result.py ---'
if [ -f /scratch/result.py ]; then sed -n '1,240p' /scratch/result.py; else echo 'NO TEMPLATE'; fi
