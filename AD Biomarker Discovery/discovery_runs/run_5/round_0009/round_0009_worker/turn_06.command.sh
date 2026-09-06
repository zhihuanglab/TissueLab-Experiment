set -euo pipefail
echo 'shared_analysis __init__.py:'
sed -n '1,220p' /shared/lib/shared_analysis/__init__.py
echo
echo 'sea_ad_lfb.py head:'
sed -n '1,260p' /shared/lib/shared_analysis/sea_ad_lfb.py
echo
echo 'stats.py head:'
sed -n '1,260p' /shared/lib/shared_analysis/stats.py
echo
echo 'zarr sample listing:'
python - <<'PY'
from pathlib import Path
p=Path('/data/H19.33.004-A12-LFB.svs.zarr')
for child in sorted(p.iterdir()):
    print(child.name, 'dir' if child.is_dir() else 'file')
    if child.is_dir():
        try:
            for grand in sorted(child.iterdir()):
                print('  ', grand.name)
        except Exception as e:
            print('   <err>', e)
PY
