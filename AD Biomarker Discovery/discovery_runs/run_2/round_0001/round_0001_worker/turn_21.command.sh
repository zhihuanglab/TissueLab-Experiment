echo '--- results.json ---'
cat /scratch/results.json
echo
echo '--- report.md ---'
cat /scratch/report.md
echo
echo '--- donor_feature_table.csv head ---'
sed -n '1,5p' /scratch/donor_feature_table.csv
echo
echo '--- replay config ---'
cat /scratch/result_replay.json
