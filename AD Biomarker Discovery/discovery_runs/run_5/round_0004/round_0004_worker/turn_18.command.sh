python - <<'PY'
import json
res=json.load(open('/scratch/results.json'))
print({k:res[k] for k in ['bootstrap_median_partial_r','bootstrap_sign_consistency','ci_lo','ci_hi']})
for r in res['ranked_variations']:
    print(r['name'], r['bootstrap_median_partial_r'], r['bootstrap_sign_consistency'], r['ci_lo'], r['ci_hi'])
PY
