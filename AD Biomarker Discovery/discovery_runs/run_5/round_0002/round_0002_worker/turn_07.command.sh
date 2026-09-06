set -euo pipefail
sed -n '1,260p' /data/autoresearch_runs/run_20260412_150614_9b299b/shared/worker_chains/ca1_pyramidal_reactive_astro_niche/sandbox/result.py
echo '--- more ---'
sed -n '261,560p' /data/autoresearch_runs/run_20260412_150614_9b299b/shared/worker_chains/ca1_pyramidal_reactive_astro_niche/sandbox/result.py
echo '--- results json ---'
cat /data/autoresearch_runs/run_20260412_150614_9b299b/shared/worker_chains/ca1_pyramidal_reactive_astro_niche/sandbox/results.json
