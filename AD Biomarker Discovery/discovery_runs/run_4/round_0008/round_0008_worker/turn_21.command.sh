set -euo pipefail
python - <<'PY'
import importlib.util, pandas as pd, numpy as np, json
spec = importlib.util.spec_from_file_location("resmod", "/scratch/result.py")
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)
with open('/scratch/context_bundle.json') as f:
    brief=json.load(f)['worker_brief']
specs=mod._variation_specs_from_brief(brief)
donor_table=mod._build_donor_table(mod.Path('/data'), brief, specs)
best=[s for s in specs if s['name']=='candidate_variant_b'][0]
loo_r,pred=mod._loo_predictions(donor_table,best['feature_column'])
pred['abs_error']=(pred['outcome']-pred['predicted']).abs()
merged=pred.merge(donor_table[['donor_id','ca1_pyramidal_n','ca1_reactive_astrocyte_n','ca1_oligodendrocyte_n',best['feature_column'],'max_age_vis','braak_numeric','cerad_ordinal','sex_binary']], on='donor_id', how='left')
print('loo_r',loo_r)
print(merged.sort_values('abs_error', ascending=False).head(10).to_string(index=False))
PY
