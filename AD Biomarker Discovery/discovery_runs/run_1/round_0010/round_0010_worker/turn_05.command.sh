set -euo pipefail
find /shared -maxdepth 3 -type d | sed -n '1,80p'
