set -euo pipefail
python /scratch/result.py > /scratch/result_stdout.txt
tail -n 120 /scratch/result_stdout.txt
