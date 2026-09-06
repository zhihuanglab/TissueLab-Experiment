set -euo pipefail
find /shared/lib -maxdepth 3 -type f | sed -n '1,120p'
