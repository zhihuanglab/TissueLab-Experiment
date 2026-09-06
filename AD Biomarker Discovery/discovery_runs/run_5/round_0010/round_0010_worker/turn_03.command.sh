set -euo pipefail
echo '--- shared_analysis.sea_ad_lfb helper signatures ---'
python - <<'PY'
import inspect
from shared_analysis import sea_ad_lfb
for name in ['build_cell_table','load_training_cohort','load_region_polygons','compute_contour_geometry']:
    fn=getattr(sea_ad_lfb,name)
    print(f'\n{name}{inspect.signature(fn)}')
    src=inspect.getsource(fn)
    print('\n'.join(src.splitlines()[:80]))
PY
echo
echo '--- previous round 4/7/9 scripts head and relevant compute functions ---'
for f in \
/Volumes/SSK\ SSD/TissueLab\ Revision/LFB/Seed48_Split/Training/autoresearch_runs/run_20260421_052101_efbee1/round_0004/round_0004_worker/result.py \
/Volumes/SSK\ SSD/TissueLab\ Revision/LFB/Seed48_Split/Training/autoresearch_runs/run_20260421_052101_efbee1/round_0007/round_0007_worker/result.py \
/Volumes/SSK\ SSD/TissueLab\ Revision/LFB/Seed48_Split/Training/autoresearch_runs/run_20260421_052101_efbee1/round_0009/round_0009_worker/result.py
do
  echo "--- $f ---"
  sed -n '1,240p' "$f" | head -n 220
done
