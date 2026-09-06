set -euo pipefail
find /data/autoresearch_runs -path '*round_0006*result.py' -o -path '*round_0006*results.json' | sed -n '1,40p'
