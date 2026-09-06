set -euo pipefail
python /scratch/result.py > /scratch/result_stdout.txt
sed -n '1,220p' /scratch/result_stdout.txt
