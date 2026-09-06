set -euo pipefail
echo '--- context_bundle.json ---'
python - <<'PY'
import json, textwrap
path='/scratch/context_bundle.json'
with open(path) as f:
    data=json.load(f)
print(json.dumps(data, indent=2)[:20000])
PY
echo
echo '--- program.md ---'
sed -n '1,220p' /scratch/program.md
echo
echo '--- result.py template ---'
sed -n '1,260p' /scratch/result.py
