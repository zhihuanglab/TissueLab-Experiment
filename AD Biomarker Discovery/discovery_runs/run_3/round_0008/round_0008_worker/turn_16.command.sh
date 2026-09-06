set -euo pipefail
grep -R "adjusted_score" -n /shared/lib /data/autoresearch_runs | head -30
