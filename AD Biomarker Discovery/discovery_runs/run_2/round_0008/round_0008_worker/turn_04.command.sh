set -euo pipefail
ls -1 /data | sed -n '1,20p'
echo
find /data/autoresearch_runs -maxdepth 4 -path '*/round_0007_worker/*' | sed -n '1,40p'
