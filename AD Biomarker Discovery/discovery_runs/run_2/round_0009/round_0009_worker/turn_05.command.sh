set -euo pipefail
echo 'shared dirs:'
find /shared -maxdepth 2 -type d | sed -n '1,80p'
echo 'shared lib files:'
find /shared/lib -maxdepth 2 -type f | sed -n '1,120p'
