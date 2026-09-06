set -e
echo '--- round 9 results.json ---'
python - <<'PY'
import json, pathlib
p=pathlib.Path('/data/autoresearch_runs/run_20260420_215745_90b105/round_0009/round_0009_worker/results.json')
obj=json.load(open(p))
print(json.dumps(obj, indent=2)[:5000])
PY
