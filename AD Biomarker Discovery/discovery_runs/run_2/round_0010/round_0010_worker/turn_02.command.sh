set -e
echo '--- worker_brief.json ---'
cat /scratch/worker_brief.json
echo
echo '--- context keys and round brief extracts ---'
python - <<'PY'
import json
p='/scratch/context_bundle.json'
with open(p) as f:
    obj=json.load(f)
print('top keys:', list(obj)[:20])
for k in ['worker_brief','candidate','round_plan','accepted_panel','current_panel','hypothesis','family']:
    if k in obj:
        print(f'\n== {k} ==')
        print(json.dumps(obj[k], indent=2)[:8000])
PY
echo
echo '--- shared embedding template head ---'
sed -n '1,260p' /shared/templates/worker_embedding_mechanistic_template.py
