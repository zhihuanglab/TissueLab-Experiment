set -euo pipefail
grep -R "Adjusted Score" -n /shared/lib /shared/templates /data/autoresearch_runs 2>/dev/null | sed -n '1,60p'
echo '---'
grep -R "penalty=" -n /shared/lib /shared/templates /data/autoresearch_runs 2>/dev/null | sed -n '1,120p'
