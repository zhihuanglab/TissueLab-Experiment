set -euo pipefail
sed -n '1,260p' /data/autoresearch_runs/run_20260421_052101_efbee1/round_0007/round_0007_worker/result.py
echo '--- results ---'
python - <<'PY'
import json
p='/data/autoresearch_runs/run_20260421_052101_efbee1/round_0007/round_0007_worker/results.json'
with open(p) as f: obj=json.load(f)
import pprint
pprint.pp(obj)
PY
