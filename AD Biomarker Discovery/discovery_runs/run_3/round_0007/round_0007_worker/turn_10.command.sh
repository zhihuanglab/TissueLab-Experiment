set -euo pipefail
grep -R "Adjusted Score" -n /data/autoresearch_runs/run_20260421_001115_585b70/round_*/*/result.py | sed -n '1,20p'
