set -e
python - <<'PY'
import sys
sys.path.insert(0,'/shared/lib')
from shared_analysis import load_training_cohort
cohort=load_training_cohort('/data')
print(cohort.columns.tolist())
print(cohort[['donor_id','cognitive_status','overall_ad_neuropath_change']].head().to_string(index=False))
PY
