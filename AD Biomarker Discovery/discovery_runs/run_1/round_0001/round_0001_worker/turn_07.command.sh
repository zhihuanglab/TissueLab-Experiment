set -euo pipefail
sed -n '1,320p' /shared/lib/shared_analysis/sea_ad_lfb.py | sed -n '1,320p'
echo '--- stats.py relevant functions ---'
grep -n "def partial_correlation\|def leave_one_out_summary\|def residualized_loo_predictive_correlation" -n /shared/lib/shared_analysis/stats.py
sed -n '1,260p' /shared/lib/shared_analysis/stats.py
