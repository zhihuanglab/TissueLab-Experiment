set -euo pipefail
python /scratch/result.py > /scratch/run_stdout.txt
sed -n '1,120p' /scratch/run_stdout.txt
echo '--- updated report ---'
sed -n '1,240p' /scratch/report.md
