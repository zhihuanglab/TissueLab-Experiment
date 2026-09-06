set -euo pipefail
python - <<'PY'
import json
from pathlib import Path
ctx=json.load(open('/scratch/context_bundle.json'))
print('TOP-LEVEL KEYS:', list(ctx))
for k in ['worker_brief','round','candidate','accepted_panel','current_panel','planned_variations','hypothesis_family','seed_candidate']:
    if k in ctx:
        print(f'\n## {k}')
        print(json.dumps(ctx[k], indent=2)[:12000])
PY

echo
echo '--- inspect one slide zarr tree (top few levels) ---'
python - <<'PY'
import zarr, os, json
path='/data/H19.33.004-A12-LFB.svs.zarr'
root=zarr.open(path, mode='r')
def walk(group, prefix='', depth=0, maxdepth=3):
    if depth>maxdepth: return
    for k,v in group.groups():
        print(prefix + k + '/')
        walk(v, prefix + '  ', depth+1, maxdepth)
    for k,v in group.arrays():
        print(prefix + f'{k}: shape={v.shape} dtype={v.dtype}')
walk(root)
PY
