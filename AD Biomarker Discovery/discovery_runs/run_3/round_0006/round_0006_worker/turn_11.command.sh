set -euo pipefail
for r in 3 4 5; do
  base="/data/autoresearch_runs/run_20260421_001115_585b70/round_000${r}/round_000${r}_worker"
  echo "--- round ${r} results.json ---"
  sed -n '1,220p' "$base/results.json"
  echo
done
