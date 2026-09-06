set -euo pipefail
grep -R "Adjusted Score" -n /shared /data /scratch 2>/dev/null | sed -n '1,80p'
