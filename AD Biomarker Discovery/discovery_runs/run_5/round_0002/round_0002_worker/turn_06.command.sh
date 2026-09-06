set -euo pipefail
grep -R "Adjusted Score" -n /shared/lib /shared/templates /data/autoresearch_runs 2>/dev/null | head -20 || true
grep -R "selection_score" -n /shared/lib /data/autoresearch_runs 2>/dev/null | head -50 || true
