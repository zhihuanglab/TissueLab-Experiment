set -euo pipefail
grep -R "adjusted_score" -n /shared/lib /shared/templates | sed -n '1,200p'
