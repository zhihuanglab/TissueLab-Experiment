set -e
echo '--- context_bundle.json ---'
python - <<'PY'
import json, textwrap
with open('/scratch/context_bundle.json') as f:
    data=json.load(f)
print(json.dumps(data, indent=2)[:30000])
PY
echo
echo '--- program.md ---'
sed -n '1,240p' /scratch/program.md
echo
echo '--- result.py template ---'
sed -n '1,260p' /scratch/result.py
