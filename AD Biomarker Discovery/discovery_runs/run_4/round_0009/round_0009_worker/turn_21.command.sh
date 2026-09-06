set -euo pipefail
nl -ba /scratch/result.py | sed -n '520,590p'
