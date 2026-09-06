set -euo pipefail
echo '--- worker_brief.json ---'
python - <<'PY'
import json
with open('/scratch/worker_brief.json') as f:
    obj=json.load(f)
print(json.dumps(obj, indent=2)[:20000])
PY
echo
echo '--- context keys + hypothesis snippet ---'
python - <<'PY'
import json, re
with open('/scratch/context_bundle.json') as f:
    obj=json.load(f)
print("top keys:", list(obj.keys()))
txt=json.dumps(obj)
for key in ['hypothesis','candidate','round','family','variation','accepted_panel','brief']:
    if key in obj:
        print(key, obj[key])
print("matches:")
for pat in ['round','candidate','hypothesis','variation','family','accepted_panel','worker_brief']:
    m=re.finditer(pat, txt, re.I)
    pos=[mm.start() for mm in m][:10]
    if pos: print(pat, pos)
PY
echo
echo '--- shared mechanistic template head ---'
sed -n '1,260p' /shared/templates/worker_embedding_mechanistic_template.py
