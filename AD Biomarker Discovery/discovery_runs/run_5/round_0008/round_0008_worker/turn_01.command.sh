set -euo pipefail
echo '--- context_bundle.json ---'
python - <<'PY'
import json, textwrap
p='/scratch/context_bundle.json'
with open(p) as f:
    obj=json.load(f)
print(json.dumps(obj, indent=2)[:30000])
PY
echo
echo '--- program.md ---'
sed -n '1,240p' /scratch/program.md || true
echo
echo '--- seeded result.py ---'
sed -n '1,260p' /scratch/result.py || true
