set -e
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

echo '--- preseed result.py head ---'
sed -n '1,220p' /scratch/result.py
