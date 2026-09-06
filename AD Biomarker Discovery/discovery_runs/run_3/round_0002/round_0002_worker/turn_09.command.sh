set -euo pipefail
sed -n '1,260p' /shared/lib/shared_analysis/stats.py
echo '---'
sed -n '1,260p' /shared/lib/shared_analysis/artifacts.py
echo '---'
sed -n '1,260p' /shared/lib/shared_analysis/sea_ad_lfb.py
