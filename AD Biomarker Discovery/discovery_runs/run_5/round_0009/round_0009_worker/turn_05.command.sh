set -euo pipefail
echo 'Shared top-level:'
find /shared -maxdepth 2 -type d | sort | sed -n '1,40p'
echo
echo 'Search for sea_ad_lfb.py and templates:'
find /shared -maxdepth 4 \( -name 'sea_ad_lfb.py' -o -name 'worker_embedding_mechanistic_template.py' -o -name 'worker_analysis_template.py' -o -name 'stats.py' \) | sort
echo
echo 'One zarr tree sample:'
find /data/H19.33.004-A12-LFB.svs.zarr -maxdepth 3 | sed -n '1,80p'
