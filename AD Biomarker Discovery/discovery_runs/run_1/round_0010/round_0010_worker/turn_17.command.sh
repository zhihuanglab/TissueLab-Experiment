set -euo pipefail
nl -ba /scratch/result.py | sed -n '326,344p'
