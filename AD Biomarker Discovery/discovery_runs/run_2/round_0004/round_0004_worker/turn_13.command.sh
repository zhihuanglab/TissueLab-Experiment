set -euo pipefail
grep -R "Adjusted Score" -n /data/autoresearch_runs | sed -n '1,20p'
