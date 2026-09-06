python - <<'PY'
import json, pandas as pd, numpy as np
from pathlib import Path
cohort=pd.read_csv('/data/training_cohort.csv')
feat=pd.read_csv('/scratch/donor_feature_table.csv')
best=json.loads(Path('/scratch/results.json').read_text())
best_col=best['feature_column']
# reconstruct LOO table by parsing stdout? better read from result? no. use run_stdout not ideal.
# Let's import result.py and compute? Use its helper.
import importlib.util, sys
spec=importlib.util.spec_from_file_location('resmod','/scratch/result.py')
mod=importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)
feature_df=feat.merge(cohort[['donor_id','sex','cognitive_status','overall_ad_neuropath_change','braak_numeric','cerad_ordinal']], on='donor_id', how='left')
loo_r, loo_table=mod._residualized_loo_table(feature_df.assign(sex_binary=feature_df['sex'].map({'Female':1.0,'Male':0.0})), best_col)
joined=loo_table.merge(feature_df[['donor_id',best_col,'sex','cognitive_status','overall_ad_neuropath_change','braak_numeric','cerad_ordinal']], on='donor_id', how='left')
joined['signed_error']=joined['predicted']-joined['actual_residualized']
print('Best col', best_col, 'loo_r', loo_r)
print('Top abs errors:')
print(joined.sort_values('abs_error', ascending=False)[['donor_id','outcome','actual_residualized','predicted','abs_error',best_col,'sex','cognitive_status','overall_ad_neuropath_change','braak_numeric','cerad_ordinal']].head(10).to_string(index=False))
print('\nFeature summary:')
print(feature_df[best_col].describe())
print('\nZero-feature donors:', int((feature_df[best_col]==0).sum()))
print(feature_df.loc[feature_df[best_col]==0,['donor_id','slope_zmem0','sex','cognitive_status','overall_ad_neuropath_change']].head(20).to_string(index=False))
PY
