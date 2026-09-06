set -euo pipefail
echo '--- prior round 6 result.py ---'
sed -n '1,260p' /data/autoresearch_runs/run_20260421_001115_585b70/round_0006/round_0006_worker/result.py
echo
echo '--- prior round 3 result.py ---'
sed -n '1,260p' /data/autoresearch_runs/run_20260421_001115_585b70/round_0003/round_0003_worker/result.py
echo
echo '--- sample zarr tree (top 3 levels) ---'
find /data/H19.33.004-A12-LFB.svs.zarr -maxdepth 3 | sed -n '1,120p'
