set -euo pipefail
export PYTHONPATH=/shared/lib:${PYTHONPATH:-}
echo '--- shared_analysis.sea_ad_lfb helper signatures ---'
python - <<'PY'
import inspect
from shared_analysis import sea_ad_lfb
for name in ['build_cell_table','load_training_cohort','load_region_polygons','compute_contour_geometry']:
    fn=getattr(sea_ad_lfb,name)
    print(f'\n{name}{inspect.signature(fn)}')
    src=inspect.getsource(fn)
    print('\n'.join(src.splitlines()[:70]))
PY
echo
echo '--- previous round 4/7/9 scripts grep relevant definitions ---'
python - <<'PY'
from pathlib import Path
files = [
Path('/Volumes/SSK SSD/TissueLab Revision/LFB/Seed48_Split/Training/autoresearch_runs/run_20260421_052101_efbee1/round_0004/round_0004_worker/result.py'),
Path('/Volumes/SSK SSD/TissueLab Revision/LFB/Seed48_Split/Training/autoresearch_runs/run_20260421_052101_efbee1/round_0007/round_0007_worker/result.py'),
Path('/Volumes/SSK SSD/TissueLab Revision/LFB/Seed48_Split/Training/autoresearch_runs/run_20260421_052101_efbee1/round_0009/round_0009_worker/result.py'),
]
for p in files:
    print(f"\n--- {p} ---")
    txt = p.read_text()
    for marker in ['FEATURE_NAME', 'def compute_donor_score', 'cKDTree', 'Reactive Astrocyte', 'contour_area']:
        idx = txt.find(marker)
        if idx >= 0:
            start = max(0, idx-300)
            end = min(len(txt), idx+2200)
            print(txt[start:end])
            print('\n[...]\n')
            break
PY
