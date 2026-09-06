set -e
grep -R "Adjusted Score" -n /shared /scratch /data/autoresearch_runs 2>/dev/null | sed -n '1,40p'
