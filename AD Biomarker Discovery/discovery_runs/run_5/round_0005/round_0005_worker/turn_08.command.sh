set -euo pipefail
find /data/autoresearch_runs -maxdepth 4 \( -name 'results.json' -o -name 'report.md' \) | sed -n '1,20p'
