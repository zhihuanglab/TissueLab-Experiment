set -euo pipefail
find /data/autoresearch_runs -path '*/round_0007_worker/sandbox/result.py' | tail -n 5
find /data/autoresearch_runs -path '*/round_0007_worker/sandbox/results.json' | tail -n 5
