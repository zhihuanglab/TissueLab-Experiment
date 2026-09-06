set -e
echo '--- round 9 result.py ---'
sed -n '1,260p' "$(python - <<'PY'
import json
obj=json.load(open('/scratch/context_bundle.json'))
print(obj['worker_brief']['accepted_panel']['members'][-1]['result_path'])
PY
)"
echo
echo '--- round 8 result.py ---'
sed -n '1,220p' "$(python - <<'PY'
import json
obj=json.load(open('/scratch/context_bundle.json'))
for m in obj['worker_brief']['accepted_panel']['members']:
    if m['round_id']==8:
        print(m['result_path'])
        break
PY
)"
echo
echo '--- shared_analysis exports and hotspot-related helpers ---'
python - <<'PY'
import shared_analysis, inspect
print('shared_analysis file:', shared_analysis.__file__)
print('exports subset:', [x for x in dir(shared_analysis) if 'region' in x.lower() or 'cell' in x.lower() or 'stat' in x.lower()][:80])
PY
python - <<'PY'
from shared_analysis import build_cell_table
import inspect
print(inspect.signature(build_cell_table))
print(inspect.getsource(build_cell_table)[:4000])
PY
