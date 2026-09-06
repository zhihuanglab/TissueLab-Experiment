set -euo pipefail
grep -R "Adjusted Score\|penalty=" -n /shared/lib /data/autoresearch_runs /scratch 2>/dev/null | sed -n '1,120p'
