set -euo pipefail
grep -R "Adjusted Score" -n /shared /scratch /data/autoresearch_runs 2>/dev/null | head -20
grep -R "is_loo_gap" -n /shared /scratch /data/autoresearch_runs 2>/dev/null | head -40
