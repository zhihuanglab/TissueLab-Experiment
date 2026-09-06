set -euo pipefail
echo '--- round4 results.json ---'
python - <<'PY'
import json
path='/data/autoresearch_runs/run_20260421_052101_efbee1/round_0004/round_0004_worker/results.json'
with open(path) as f: obj=json.load(f)
for k,v in obj.items():
    if k not in ('ranked_variations','per_donor_loo'):
        print(k,':',v)
print('ranked top:', obj.get('ranked_variations',[])[:3])
print('per_donor_loo sample:', obj.get('per_donor_loo',[])[:3])
PY
echo '--- round4 report ---'
sed -n '1,220p' /data/autoresearch_runs/run_20260421_052101_efbee1/round_0004/round_0004_worker/report.md
echo '--- round4 result.py ---'
sed -n '1,260p' /data/autoresearch_runs/run_20260421_052101_efbee1/round_0004/round_0004_worker/result.py
