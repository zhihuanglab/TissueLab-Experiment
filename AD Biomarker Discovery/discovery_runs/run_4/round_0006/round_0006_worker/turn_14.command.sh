set -euo pipefail
echo '--- round5 results.json ---'
python - <<'PY'
import json
path='/data/autoresearch_runs/run_20260421_014751_94b13b/round_0005/round_0005_worker/results.json'
with open(path) as f:
    obj=json.load(f)
print(json.dumps(obj, indent=2)[:4000])
PY
echo
echo '--- round5 result.py head ---'
sed -n '1,260p' /data/autoresearch_runs/run_20260421_014751_94b13b/round_0005/round_0005_worker/result.py
