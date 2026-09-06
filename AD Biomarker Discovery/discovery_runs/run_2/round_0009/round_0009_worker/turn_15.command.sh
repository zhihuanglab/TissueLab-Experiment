set -euo pipefail
ls -l /scratch/result.py /scratch/results.json /scratch/report.md /scratch/donor_feature_table.csv /scratch/feature_spec.json
echo '--- results.json summary ---'
python - <<'PY'
import json
with open('/scratch/results.json') as f:
    r=json.load(f)
print({k:r.get(k) for k in ['best_variation','feature_column','partial_r','selection_score','loo_predictive_r','is_loo_gap','gap_penalty','adjusted_score']})
print('ranked_variations len', len(r.get('ranked_variations', [])))
print('best feature column present in table?', r['feature_column'] in open('/scratch/donor_feature_table.csv').readline().strip().split(','))
PY
echo '--- report.md ---'
sed -n '1,220p' /scratch/report.md
