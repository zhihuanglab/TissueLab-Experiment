set -e
python - <<'PY'
import json, pandas as pd
from pathlib import Path
res = json.loads(Path('/scratch/results.json').read_text())
loo = pd.DataFrame(res['best_variation_loo']).sort_values('absolute_error', ascending=False)
cohort = pd.read_csv('/data/training_cohort.csv')
top = loo.head(8).merge(cohort, on=['donor_id','slide_name'], how='left')
print(top[['donor_id','outcome','predicted','absolute_error','ca1_pyramidal_with_reactive_astro_within_80px_fraction','sex','max_age_vis','braak_numeric','cerad_ordinal']].to_string(index=False))
PY
