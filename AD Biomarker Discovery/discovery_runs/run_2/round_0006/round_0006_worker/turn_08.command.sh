set -euo pipefail
echo '--- results.json ---'
sed -n '1,220p' /scratch/results.json
echo '--- report.md ---'
sed -n '1,240p' /scratch/report.md
