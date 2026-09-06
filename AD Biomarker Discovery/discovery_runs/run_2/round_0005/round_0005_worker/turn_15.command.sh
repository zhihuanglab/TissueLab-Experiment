set -euo pipefail
python /scratch/result.py | tee /scratch/run_stdout.txt
ls -l /scratch/results.json /scratch/report.md /scratch/donor_feature_table.csv
