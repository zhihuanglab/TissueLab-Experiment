set -euo pipefail
echo '--- worker_brief.json ---'
cat /scratch/worker_brief.json
echo
echo '--- seeded result.py ---'
sed -n '1,240p' /scratch/result.py
