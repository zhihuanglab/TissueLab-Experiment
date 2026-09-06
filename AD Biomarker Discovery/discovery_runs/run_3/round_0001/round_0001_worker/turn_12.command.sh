set -euo pipefail
grep -n "CANONICAL_VARIATION_NAME\|CANONICAL_FEATURE_COLUMN" /scratch/result.py | sed -n '1,20p'
