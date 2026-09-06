set -euo pipefail
echo '--- context_bundle.json ---'
python - <<'PY'
import json, textwrap
with open('/scratch/context_bundle.json') as f:
    obj=json.load(f)
print(json.dumps(obj, indent=2)[:20000])
PY
echo
echo '--- program.md ---'
sed -n '1,240p' /scratch/program.md
echo
echo '--- result.py template head ---'
sed -n '1,240p' /scratch/result.py
