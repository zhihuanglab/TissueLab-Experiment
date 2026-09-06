set -euo pipefail
echo '--- context_bundle.json ---'
cat /scratch/context_bundle.json
echo
echo '--- program.md ---'
cat /scratch/program.md
echo
echo '--- result.py template head ---'
sed -n '1,220p' /scratch/result.py
