set -euo pipefail
grep -n "def write_donor_feature_table\|def write_results_payload" -n /shared/lib/shared_analysis/artifacts.py
sed -n '260,380p' /shared/lib/shared_analysis/artifacts.py
