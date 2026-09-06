set -euo pipefail
find /data/autoresearch_runs -maxdepth 4 \( -name results.json -o -name report.md -o -name result.py \) | sed -n '1,60p'
