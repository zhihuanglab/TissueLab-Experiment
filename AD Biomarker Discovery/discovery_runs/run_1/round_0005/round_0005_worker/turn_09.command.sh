set -e
echo 'search for Adjusted Score / penalty'
grep -R "Adjusted Score" -n /shared/lib /shared/templates /data/autoresearch_runs 2>/dev/null | sed -n '1,40p' || true
grep -R "is_loo_gap" -n /shared/lib /shared/templates /data/autoresearch_runs 2>/dev/null | sed -n '1,80p' || true
