set -euo pipefail
echo '--- prior round 6 result.py ---'
sed -n '1,260p' "/data/autoresearch_runs/run_20260406_110805_82cf2b/round_0006/Worker 1/candidate_bundle/artifacts/result.py" || true
echo
echo '--- prior round 6 results.json ---'
python - <<'PY'
import json, glob, os
paths=glob.glob('/data/autoresearch_runs/run_20260406_110805_82cf2b/round_0006/**/results.json', recursive=True)
print(paths[:5])
for p in paths[:3]:
    print('PATH',p)
    with open(p) as f:
        obj=json.load(f)
    print(json.dumps(obj, indent=2)[:4000])
    print()
PY
