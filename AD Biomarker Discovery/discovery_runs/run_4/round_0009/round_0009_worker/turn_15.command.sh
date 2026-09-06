set -euo pipefail
python - <<'PY'
import pandas as pd
co=pd.read_csv('/data/training_cohort.csv')
print(co.columns.tolist())
print(co[['donor_id','cognitive_status','overall_ad_neuropath_change']].head(10).to_string())
print(co['cognitive_status'].value_counts(dropna=False))
print(co['overall_ad_neuropath_change'].value_counts(dropna=False))
PY
