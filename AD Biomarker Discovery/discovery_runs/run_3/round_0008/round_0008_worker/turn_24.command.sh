set -euo pipefail
python - <<'PY'
import importlib.util, json
spec=importlib.util.spec_from_file_location('result','/scratch/result.py')
mod=importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)
print('FEATURE_NAME', mod.FEATURE_NAME)
print('FEATURE_COLUMN', mod.FEATURE_COLUMN)
print('score H20.33.018', mod.compute_donor_score(donor_id='H20.33.018', data_root='/data'))
PY
