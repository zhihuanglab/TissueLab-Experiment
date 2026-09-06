set -euo pipefail
python - <<'PY'
import scipy, tifffile, numpy
print('ok', scipy.__version__, tifffile.__version__)
PY
