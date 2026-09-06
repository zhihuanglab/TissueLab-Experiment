set -euo pipefail
grep -n "def build_cell_table\|center\|area =" -n /shared/lib/shared_analysis/sea_ad_lfb.py | sed -n '1,220p'
