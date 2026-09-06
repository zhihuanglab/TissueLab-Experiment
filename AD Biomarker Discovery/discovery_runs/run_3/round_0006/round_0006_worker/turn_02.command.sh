set -euo pipefail
echo '--- worker_brief.json ---'
python - <<'PY'
import json
with open('/scratch/worker_brief.json') as f:
    data=json.load(f)
print(json.dumps(data, indent=2))
PY
echo
echo '--- shared embedding template head ---'
sed -n '1,260p' /shared/templates/worker_embedding_mechanistic_template.py
