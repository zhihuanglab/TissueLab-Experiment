set -euo pipefail
echo '--- patched constants ---'
grep -n 'CANONICAL_VARIATION\|FEATURE_NAME\|FEATURE_COLUMN' /scratch/result.py | sed -n '1,10p'
echo '--- results.json ---'
cat /scratch/results.json
echo
echo '--- report.md ---'
sed -n '1,220p' /scratch/report.md
