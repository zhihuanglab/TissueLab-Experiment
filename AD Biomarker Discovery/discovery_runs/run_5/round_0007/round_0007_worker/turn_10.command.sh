set -euo pipefail
grep -R "Adjusted Score" -n /shared /data/autoresearch_runs 2>/dev/null | head -20
grep -R "\"adjusted_score\"" -n /shared /data/autoresearch_runs 2>/dev/null | head -20
