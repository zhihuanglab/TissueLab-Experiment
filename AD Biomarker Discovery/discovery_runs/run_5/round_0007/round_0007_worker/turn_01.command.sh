set -euo pipefail
echo '--- context_bundle.json ---'
python - <<'PY'
import json, pathlib
p=pathlib.Path('/scratch/context_bundle.json')
obj=json.loads(p.read_text())
import pprint
pprint.pp(obj)
PY
echo '--- program.md ---'
sed -n '1,220p' /scratch/program.md
echo '--- existing result.py template ---'
sed -n '1,260p' /scratch/result.py || true
