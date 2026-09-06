set -euo pipefail
for r in 4 5 6; do
  p=$(find /data/autoresearch_runs -path "*/round_$(printf '%04d' $r)/*/result.py" | head -1 || true)
  echo "--- round $r: $p ---"
  [ -n "$p" ] && sed -n '1,260p' "$p" | sed -n '1,220p'
done
