set -euo pipefail
python - <<'PY'
import sys, inspect, textwrap
sys.path.insert(0, '/shared/lib')
import shared_analysis.artifacts as art
print(textwrap.dedent(inspect.getsource(art.resolve_covariate_names))[:2000])
print(textwrap.dedent(inspect.getsource(art.validate_donor_feature_table_columns))[:2000])
print('required fields', art.RESULTS_REQUIRED_FIELDS)
PY
