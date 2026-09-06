set -euo pipefail
grep -n "RESULTS_REQUIRED_FIELDS" -n /shared/lib/shared_analysis/artifacts.py
sed -n '1,120p' /shared/lib/shared_analysis/artifacts.py
