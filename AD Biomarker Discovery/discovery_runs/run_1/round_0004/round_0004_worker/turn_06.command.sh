set -euo pipefail
find /shared -maxdepth 4 \( -type d -name 'shared_analysis' -o -type f -name 'shared_analysis.py' -o -type f -name '*.py' \) | sed -n '1,200p'
