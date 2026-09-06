set -euo pipefail
echo 'shared tree'
find /shared -maxdepth 3 -type d | sed -n '1,120p'
echo 'python path candidates'
find /shared -maxdepth 4 -type f \( -name '__init__.py' -o -name 'sea_ad_lfb.py' -o -name 'stats.py' \) | sed -n '1,120p'
echo 'data autoresearch runs'
find /data/autoresearch_runs -maxdepth 4 -type f \( -name 'result.py' -o -name 'results.json' -o -name 'report.md' \) | sed -n '1,80p'
