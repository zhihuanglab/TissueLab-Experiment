set -euo pipefail
find /data/autoresearch_runs -maxdepth 4 -path '*/round_0007/*' -type f | sed -n '1,40p'
