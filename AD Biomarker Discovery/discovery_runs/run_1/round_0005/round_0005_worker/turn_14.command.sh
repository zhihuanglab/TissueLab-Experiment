set -e
python - <<'PY'
import sys
sys.path.insert(0,'/shared/lib')
import shared_analysis
print('DEFAULT_CONFOUNDS', shared_analysis.DEFAULT_CONFOUNDS)
PY
