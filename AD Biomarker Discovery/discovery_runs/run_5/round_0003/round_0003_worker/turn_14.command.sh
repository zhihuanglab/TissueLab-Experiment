python - <<'PY'
import json, pandas as pd, importlib.util
from pathlib import Path
cohort=pd.read_csv('/data/training_cohort.csv')
feat=pd.read_csv('/scratch/donor_feature_table.csv')
best=json.loads(Path('/scratch/results.json').read_text())
best_col=best['feature_column']

spec=importlib.util.spec_from_file_location('resmod','/scratch/result.py')
mod=importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

df=feat.copy()
loo_r, loo_table=mod._residualized_loo_table(df, best_col)
joined=loo_table.merge(cohort[['donor_id','max_age_vis','braak_numeric','cerad_ordinal','sex']], on='donor_id', how='left')
print('Top abs errors:')
print(joined.sort_values('abs_error', ascending=False)[['donor_id','outcome','actual_residualized','predicted','abs_error',best_col,'max_age_vis','braak_numeric','cerad_ordinal','sex']].head(12).to_string(index=False))
print('\nFeature quantiles:')
print(df[best_col].quantile([0,0.1,0.25,0.5,0.75,0.9,1.0]).to_string())
print('\nHighest feature donors:')
print(df[['donor_id','slope_zmem0',best_col]].sort_values(best_col, ascending=False).head(10).to_string(index=False))
print('\nZero-feature donors count', int((df[best_col]==0).sum()))
print(df.loc[df[best_col]==0,['donor_id','slope_zmem0',best_col,'max_age_vis','braak_numeric','cerad_ordinal','sex']].to_string(index=False))
PY
