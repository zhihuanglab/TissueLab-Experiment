set -euo pipefail
python - <<'PY'
import sys
sys.path.insert(0, '/shared/lib')
import shared_analysis.artifacts as art
print(art._COVARIATE_ALIASES)
PY
