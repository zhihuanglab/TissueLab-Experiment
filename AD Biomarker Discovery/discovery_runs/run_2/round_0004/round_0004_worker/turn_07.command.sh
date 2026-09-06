set -euo pipefail
find /shared/lib -maxdepth 2 -type f | sed -n '1,50p'
