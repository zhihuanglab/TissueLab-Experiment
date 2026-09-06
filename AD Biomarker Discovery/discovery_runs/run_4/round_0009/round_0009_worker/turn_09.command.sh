set -euo pipefail
grep -n "def load_centroids\|def load_class_ids\|def load_class_lookup\|def load_region_polygons" -n /shared/lib/shared_analysis/sea_ad_lfb.py
sed -n '1,260p' /shared/lib/shared_analysis/sea_ad_lfb.py | sed -n '1,240p'
