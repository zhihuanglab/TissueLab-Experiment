set -euo pipefail
echo '--- /scratch/context_bundle.json ---'
python - <<'PY'
import json, pprint
with open('/scratch/context_bundle.json') as f:
    obj=json.load(f)
pprint.pp(obj)
PY
echo '--- /scratch/program.md ---'
cat /scratch/program.md
echo '--- /scratch listing ---'
ls -lah /scratch
echo '--- preseeded /scratch/result.py ---'
sed -n '1,220p' /scratch/result.py
