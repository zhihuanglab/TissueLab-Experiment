set -euo pipefail
echo '--- results.json head ---'
sed -n '1,220p' /scratch/results.json
echo
echo '--- report.md ---'
sed -n '1,260p' /scratch/report.md
