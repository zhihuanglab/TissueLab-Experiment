set -euo pipefail
echo '--- context_bundle.json ---'
python - <<'PY'
import json, textwrap
path='/scratch/context_bundle.json'
with open(path) as f:
    obj=json.load(f)
print(json.dumps(obj, indent=2)[:30000])
PY
echo
echo '--- program.md ---'
sed -n '1,240p' /scratch/program.md
echo
echo '--- existing result.py template head ---'
sed -n '1,240p' /scratch/result.py || true
