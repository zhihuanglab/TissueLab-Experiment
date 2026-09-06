set -euo pipefail
grep -n "write_donor_feature_table" -n /shared/lib/shared_analysis/artifacts.py
sed -n '260,360p' /shared/lib/shared_analysis/artifacts.py
