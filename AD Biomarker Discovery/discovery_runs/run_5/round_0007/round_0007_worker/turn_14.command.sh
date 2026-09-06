python - <<'PY'
from pathlib import Path
p = Path('/scratch/result.py')
text = p.read_text()
old = '''    analyzable = df[["donor_id", feature_col, OUTCOME_COLUMN, *CONFOUND_COLUMNS]].dropna().copy()
'''
new = '''    analyzable = df.dropna(subset=["donor_id", feature_col, OUTCOME_COLUMN, *CONFOUND_COLUMNS]).copy()
'''
if old not in text:
    raise SystemExit('target snippet not found')
p.write_text(text.replace(old, new))
PY
python /scratch/result.py > /scratch/run.log 2>&1
status=$?
echo "script exit status: $status"
sed -n '1,260p' /scratch/run.log
exit $status
