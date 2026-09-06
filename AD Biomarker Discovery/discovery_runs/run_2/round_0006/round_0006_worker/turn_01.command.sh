set -euo pipefail
echo '--- context_bundle.json ---'
python - <<'PY'
import json, textwrap
with open('/scratch/context_bundle.json') as f:
    obj=json.load(f)
print(json.dumps(obj, indent=2)[:30000])
PY
echo '--- program.md ---'
sed -n '1,220p' /scratch/program.md || true
echo '--- seeded result.py ---'
sed -n '1,260p' /scratch/result.py || true
