set -euo pipefail
echo '--- shared_analysis/__init__.py ---'
sed -n '1,220p' /shared/lib/shared_analysis/__init__.py
echo '--- sea_ad_lfb.py relevant bits ---'
grep -n "def build_cell_table\|def load_training_cohort\|def load_region_polygons\|def compute_contour_geometry" -n /shared/lib/shared_analysis/sea_ad_lfb.py
sed -n '1,260p' /shared/lib/shared_analysis/sea_ad_lfb.py
echo '--- round 5 result.py ---'
sed -n '1,260p' /data/autoresearch_runs/run_20260420_215745_90b105/round_0005/round_0005_worker/result.py
echo '--- round 3 result.py ---'
sed -n '1,260p' /data/autoresearch_runs/run_20260420_215745_90b105/round_0003/round_0003_worker/result.py
