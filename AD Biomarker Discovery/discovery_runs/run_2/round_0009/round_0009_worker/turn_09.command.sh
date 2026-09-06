set -euo pipefail
grep -n "def write_donor_feature_table" -n /shared/lib/shared_analysis/artifacts.py
sed -n '300,380p' /shared/lib/shared_analysis/artifacts.py
