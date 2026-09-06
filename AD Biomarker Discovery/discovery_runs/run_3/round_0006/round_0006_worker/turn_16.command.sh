set -euo pipefail
for r in 4 5; do
  base="/data/autoresearch_runs/run_20260421_001115_585b70/round_000${r}/round_000${r}_worker/result.py"
  echo "--- round ${r} script head ---"
  sed -n '1,260p' "$base" | sed -n '1,220p'
  echo
done
