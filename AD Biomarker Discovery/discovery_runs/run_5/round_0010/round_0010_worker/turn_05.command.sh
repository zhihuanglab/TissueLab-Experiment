set -euo pipefail
echo '--- /shared ---'
find /shared -maxdepth 3 -type d | sed -n '1,80p'
echo '--- python sys.path and candidates ---'
python - <<'PY'
import sys, os
print(sys.path)
for p in ['/shared/lib','/shared','/shared/cache','/shared/templates']:
    print(p, os.path.exists(p), os.listdir(p)[:10] if os.path.isdir(p) else None)
PY
echo '--- inspect one zarr store and csv head ---'
ls -1 /data | sed -n '1,12p'
python - <<'PY'
import pandas as pd
cohort=pd.read_csv('/data/training_cohort.csv')
print(cohort[['donor_id','slide_name','slope_zmem0','max_age_vis','braak_numeric','cerad_ordinal','sex']].head())
PY
echo '--- zarr tree head ---'
python - <<'PY'
import zarr
from pathlib import Path
p=next(Path('/data').glob('*.svs.zarr'))
print('sample zarr:', p.name)
root=zarr.open(str(p), mode='r')
print(root.tree())
PY
