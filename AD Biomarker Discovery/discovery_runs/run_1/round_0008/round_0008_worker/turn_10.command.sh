set -euo pipefail
grep -R "adjusted_score" -n /shared/lib /shared/templates /data/autoresearch_runs | head -n 30
