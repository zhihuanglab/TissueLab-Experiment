set -euo pipefail
python - <<'PY'
from pathlib import Path
paths=[
'/data/autoresearch_runs',
'/Volumes/SSK SSD/TissueLab Revision/LFB/Seed48_Split/Training/autoresearch_runs'
]
for p in paths:
    print(p, Path(p).exists())
PY
# try to inspect shared previous result if accessible
for p in \
"/Volumes/SSK SSD/TissueLab Revision/LFB/Seed48_Split/Training/autoresearch_runs/run_20260421_014751_94b13b/round_0007/round_0007_worker/result.py" \
"/Volumes/SSK SSD/TissueLab Revision/LFB/Seed48_Split/Training/autoresearch_runs/run_20260421_014751_94b13b/round_0007/round_0007_worker/results.json"
do
  echo "--- $p ---"
  sed -n '1,220p' "$p" | head -n 220
done
