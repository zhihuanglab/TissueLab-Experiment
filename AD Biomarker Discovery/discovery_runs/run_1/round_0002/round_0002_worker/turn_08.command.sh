set -euo pipefail
python - <<'PY'
import pandas as pd
df=pd.read_csv('/data/training_cohort.csv')
print(df.columns.tolist())
print(df.head(1).to_dict(orient='records')[0])
PY
# patch result.py for missing cognitive_status
python - <<'PY'
from pathlib import Path
p=Path('/scratch/result.py')
text=p.read_text()
text=text.replace('            "cognitive_status": row.cognitive_status,\n            "overall_ad_neuropath_change": row.overall_ad_neuropath_change,\n',
                  '            "cognitive_status": getattr(row, "cognitive_status", None),\n            "overall_ad_neuropath_change": getattr(row, "overall_ad_neuropath_change", None),\n')
text=text.replace('        donor_bits.append(\n            f"{row.donor_id} (outcome={row.slope_zmem0:.3f}, pred={row.predicted:.3f}, "\n            f"{best[\'feature_column\']}={getattr(row, best[\'feature_column\']):.4f}, "\n            f"status={row.cognitive_status}, AD_change={row.overall_ad_neuropath_change})"\n        )\n',
                  '        donor_bits.append(\n            f"{row.donor_id} (outcome={row.slope_zmem0:.3f}, pred={row.predicted:.3f}, "\n            f"{best[\'feature_column\']}={getattr(row, best[\'feature_column\']):.4f}, "\n            f"status={getattr(row, \'cognitive_status\', None)}, AD_change={getattr(row, \'overall_ad_neuropath_change\', None)})"\n        )\n')
p.write_text(text)
PY
python /scratch/result.py | tee /scratch/run_output.txt
echo '--- result.py canonical variation line ---'
grep 'CANONICAL_VARIATION' /scratch/result.py
echo '--- results.json ---'
python - <<'PY'
import json
obj=json.load(open('/scratch/results.json'))
print(json.dumps({k: obj[k] for k in ['best_variation','feature_column','partial_r','selection_score','loo_predictive_r','is_loo_gap','adjusted_score']}, indent=2))
PY
echo '--- report head ---'
sed -n '1,120p' /scratch/report.md
