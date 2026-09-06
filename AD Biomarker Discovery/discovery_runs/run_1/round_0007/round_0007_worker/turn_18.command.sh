set -euo pipefail
time python /scratch/result.py > /scratch/result_stdout.txt
tail -n 80 /scratch/result_stdout.txt
