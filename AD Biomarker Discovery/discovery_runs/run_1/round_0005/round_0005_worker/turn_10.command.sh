set -e
grep -R "is_loo_gap" -n /shared/lib /shared/templates 2>/dev/null || true
grep -R "Adjusted Score" -n /shared/lib /shared/templates 2>/dev/null || true
