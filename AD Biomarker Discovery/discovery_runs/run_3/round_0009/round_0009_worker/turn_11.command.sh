set -euo pipefail
grep -R "Adjusted Score" -n /shared/lib /shared/templates /data/autoresearch_runs 2>/dev/null | head -20
