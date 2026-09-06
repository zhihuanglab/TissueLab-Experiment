set -euo pipefail
python /scratch/result.py > /scratch/result_stdout.txt
sed -n '1,220p' /scratch/result_stdout.txt
python - <<'PY'
import json
d=json.load(open('/scratch/results.json'))
print('\nRESULTS KEYS', sorted(d.keys()))
print('BEST', d['best_variation'], d['feature_column'])
print('RANKED', d['ranked_variations'])
PY
