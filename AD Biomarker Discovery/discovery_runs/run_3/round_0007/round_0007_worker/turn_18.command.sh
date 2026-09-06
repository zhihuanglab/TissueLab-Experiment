set -euo pipefail
echo '--- rewritten canonical constants ---'
grep -n "AUTO_" /scratch/result.py
echo
echo '--- report.md ---'
sed -n '1,240p' /scratch/report.md
echo
echo '--- results.json head ---'
sed -n '1,220p' /scratch/results.json
