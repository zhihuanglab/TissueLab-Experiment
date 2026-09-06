set -e
echo '--- context_bundle.json ---'
cat /scratch/context_bundle.json
echo
echo '--- program.md ---'
cat /scratch/program.md
echo
echo '--- preseeded result.py ---'
if [ -f /scratch/result.py ]; then sed -n '1,220p' /scratch/result.py; else echo 'NO TEMPLATE'; fi
