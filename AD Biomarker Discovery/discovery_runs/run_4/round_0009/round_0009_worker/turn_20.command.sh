set -euo pipefail
nl -ba /scratch/result.py | sed -n '460,520p'
