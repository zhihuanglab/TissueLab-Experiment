set -euo pipefail
find /shared -maxdepth 3 \( -type d -o -type f \) | grep '/shared/lib' | sed -n '1,120p'
