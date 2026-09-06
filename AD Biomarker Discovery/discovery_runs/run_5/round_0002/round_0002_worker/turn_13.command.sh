set -euo pipefail
python /scratch/result.py | tee /scratch/run_stdout.txt
echo '--- results.json ---'
cat /scratch/results.json
echo '--- report.md ---'
cat /scratch/report.md
