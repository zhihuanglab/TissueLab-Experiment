set -euo pipefail
sed -n '1,220p' /shared/lib/shared_analysis/__init__.py
echo '--- sea_ad_lfb.py ---'
grep -n "def build_cell_table\|def load_training_cohort\|def compute_contour_geometry" -n /shared/lib/shared_analysis/sea_ad_lfb.py
sed -n '1,260p' /shared/lib/shared_analysis/sea_ad_lfb.py
echo '--- stats.py ---'
sed -n '1,260p' /shared/lib/shared_analysis/stats.py
