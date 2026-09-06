set -euo pipefail
grep -R "Adjusted Score\|penalty\|selection_score" -n /shared/lib/shared_analysis | sed -n '1,200p'
