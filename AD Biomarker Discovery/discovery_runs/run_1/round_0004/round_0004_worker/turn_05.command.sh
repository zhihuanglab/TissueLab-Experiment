set -euo pipefail
echo 'PYTHONPATH:' "${PYTHONPATH:-}"
find /shared -maxdepth 3 -type d -name 'shared_analysis' -o -maxdepth 3 -type f | sed -n '1,80p'
