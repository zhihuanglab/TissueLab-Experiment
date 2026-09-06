python - <<'PY'
import json
res=json.load(open('/scratch/results.json'))
print(json.dumps(res['ranked_variations'], indent=2))
PY
