set -euo pipefail
grep -n "sorted(metrics" -n /data/autoresearch_runs/run_20260421_001115_585b70/round_0006/round_0006_worker/result.py
grep -n "selection_score" /data/autoresearch_runs/run_20260421_001115_585b70/round_0006/round_0006_worker/result.py | sed -n '1,20p'
