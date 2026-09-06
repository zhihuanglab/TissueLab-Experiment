set -euo pipefail
grep -R "Adjusted Score" -n /data/autoresearch_runs | head -n 20
