set -euo pipefail
echo '--- report.md ---'
sed -n '1,240p' /scratch/report.md
echo
echo '--- results.json ---'
cat /scratch/results.json
