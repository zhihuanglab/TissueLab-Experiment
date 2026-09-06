set -euo pipefail
grep -R "Adjusted Score" -n /shared /data/autoresearch_runs | sed -n '1,40p' || true
grep -R "\"adjusted_score\"" -n /data/autoresearch_runs | sed -n '1,40p' || true
