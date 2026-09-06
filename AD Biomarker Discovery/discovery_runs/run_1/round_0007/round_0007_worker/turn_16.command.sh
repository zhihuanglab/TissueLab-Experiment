set -euo pipefail
grep -n "def write_results_payload\|def build_results_payload" -n /shared/lib/shared_analysis/artifacts.py
sed -n '320,420p' /shared/lib/shared_analysis/artifacts.py
