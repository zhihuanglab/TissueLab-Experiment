set -e
python - <<'PY'
import sys, json, pandas as pd
sys.path.insert(0,'/shared/lib')
from shared_analysis.stats import residualized_loo_predictive_correlation
cohort = pd.read_csv('/scratch/donor_feature_table.csv')
r = residualized_loo_predictive_correlation(cohort, feature_col='ca1_pyramidal_with_reactive_astro_within_80px_fraction', outcome_col='slope_zmem0', confounds=['max_age_vis','braak_numeric','cerad_ordinal','sex_binary'])
print(r)
PY
