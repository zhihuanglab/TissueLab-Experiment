set -euo pipefail
echo '--- worker_brief.json ---'
python - <<'PY'
import json
with open('/scratch/worker_brief.json') as f:
    obj=json.load(f)
print(json.dumps(obj, indent=2))
PY
echo
echo '--- list scratch ---'
ls -la /scratch | sed -n '1,120p'
