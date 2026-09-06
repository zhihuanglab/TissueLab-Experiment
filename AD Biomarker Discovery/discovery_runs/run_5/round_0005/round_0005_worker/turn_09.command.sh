set -euo pipefail
python - <<'PY'
import json, glob, os
matches=[]
for path in glob.glob('/data/autoresearch_runs/**/results.json', recursive=True):
    try:
        with open(path) as f: obj=json.load(f)
    except Exception:
        continue
    if obj.get('feature_name') in {
        'ca1_reactive_astro_fraction',
        'ca1_pyramidal_near_reactive_astro_fraction_30um',
        'ca1_pyramidal_reactive_astro_crowding_fraction_k3_30um',
        'ca1_peripyramidal_reactive_astro_area_median_30um_um2',
    }:
        matches.append((obj.get('feature_name'), path))
for m in matches[:20]:
    print(m)
PY
