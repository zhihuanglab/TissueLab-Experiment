set -euo pipefail
python - <<'PY'
import pandas as pd, json
res=json.load(open('/scratch/results.json'))
best=res['feature_column']
per=pd.read_csv('/scratch/donor_feature_table.csv')[['donor_id',best,'reactive_niche_pyramidal_n','ca1_corpora_amylacea_n']]
loo_rows=res['ranked_variations'][0]  # not per donor
# parse from stdout unavailable, recompute by reading report? easier load from result? no per_donor absent results
# use donor_feature_table and estimates from report? Let's read per-donor from printed? skip
PY
