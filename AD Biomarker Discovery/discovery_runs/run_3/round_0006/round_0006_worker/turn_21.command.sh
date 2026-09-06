set -euo pipefail
echo '--- results.json ---'
sed -n '1,240p' /scratch/results.json
echo
echo '--- report.md ---'
sed -n '1,260p' /scratch/report.md
