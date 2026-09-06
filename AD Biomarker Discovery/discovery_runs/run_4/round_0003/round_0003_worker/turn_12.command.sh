set -e
echo '--- canonical variation line ---'
grep -n 'CANONICAL_VARIATION' /scratch/result.py | sed -n '1,3p'
echo '--- results.json head ---'
sed -n '1,220p' /scratch/results.json
