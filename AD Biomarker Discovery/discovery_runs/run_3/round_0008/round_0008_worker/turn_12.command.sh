set -euo pipefail
find /data/autoresearch_runs -name results.json | head -20 | sed -n '1,20p'
