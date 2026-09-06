set -euo pipefail
grep -n "def partial_correlation\|def bootstrap_partial_correlation\|def leave_one_out_summary\|def residualized_loo_predictive_correlation\|def analyze_loo_instability" -n /shared/lib/shared_analysis/stats.py
sed -n '1,220p' /shared/lib/shared_analysis/stats.py
