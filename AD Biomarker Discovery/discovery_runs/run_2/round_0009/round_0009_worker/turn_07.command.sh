set -euo pipefail
grep -n "def build_cell_table" -n /shared/lib/shared_analysis/sea_ad_lfb.py
sed -n '280,430p' /shared/lib/shared_analysis/sea_ad_lfb.py
echo '--- assign helpers ---'
grep -n "def assign_centroids_to_regions" -n /shared/lib/shared_analysis/sea_ad_lfb.py
sed -n '220,300p' /shared/lib/shared_analysis/sea_ad_lfb.py
