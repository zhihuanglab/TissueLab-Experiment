set -euo pipefail
python - <<'PY'
import sys, math, json, warnings
sys.path.insert(0,'/shared/lib')
warnings.filterwarnings('ignore')
from pathlib import Path
import numpy as np, pandas as pd
from shared_analysis import load_training_cohort, build_cell_table, partial_correlation, residualized_loo_predictive_correlation

cohort = load_training_cohort('/data').copy()
cohort['sex_binary'] = cohort['sex'].map(lambda s: 1.0 if str(s).lower().startswith('m') else 0.0)
rows=[]
for _,r in cohort.iterrows():
    slide = Path('/data')/r['slide_name']
    cells = build_cell_table(slide, include_regions=True, include_geometry=False)
    valid = ~cells['cell_type'].astype(str).str.lower().isin(['unknown','unlabeled','unlabelled',''])
    out={'donor_id':r['donor_id']}
    for reg in ['CA1','CA2']:
        m = (cells['region']==reg) & valid
        out[f'{reg.lower()}_valid'] = int(m.sum())
        out[f'{reg.lower()}_pyr'] = int(((cells['cell_type']=='Pyramidal Neuron') & m).sum())
    rows.append(out)
f = pd.DataFrame(rows)
df = cohort.merge(f,on='donor_id')
df['candidate_variant_a'] = np.where(df['ca1_valid']>=100, df['ca1_pyr']/df['ca1_valid'], np.nan)
den = df['ca1_valid']+df['ca2_valid']
num = df['ca1_pyr']+df['ca2_pyr']
df['candidate_variant_b'] = np.where(den>=100, num/den, np.nan)
for col in ['candidate_variant_a','candidate_variant_b']:
    pr = partial_correlation(df, feature_col=col, outcome_col='slope_zmem0', confounds=['max_age_vis','braak_numeric','cerad_ordinal','sex_binary'])
    loo = residualized_loo_predictive_correlation(df, feature_col=col, outcome_col='slope_zmem0', confounds=['max_age_vis','braak_numeric','cerad_ordinal','sex_binary'])
    n = int(pr['n']); sel = abs(pr['partial_r'])*(n/len(df))
    gap = abs(pr['partial_r'])-abs(loo) if np.isfinite(loo) else float('nan')
    print(col, pr, 'selection', sel, 'loo', loo, 'gap', gap)
PY
