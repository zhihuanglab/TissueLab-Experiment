set -euo pipefail
echo $PYTHONPATH
find /shared -maxdepth 3 -type d | sed -n '1,60p'
find /shared -maxdepth 3 -type f | sed -n '1,80p'
