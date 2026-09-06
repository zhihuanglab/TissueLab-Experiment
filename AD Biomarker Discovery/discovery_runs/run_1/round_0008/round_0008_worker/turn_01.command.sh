set -euo pipefail
echo '--- context_bundle.json ---'
python - <<'PY'
import json
p='/scratch/context_bundle.json'
with open(p) as f:
    obj=json.load(f)
import pprint
pprint.pp(obj)
PY
echo '--- program.md ---'
cat /scratch/program.md
echo '--- seeded result.py ---'
sed -n '1,240p' /scratch/result.py
