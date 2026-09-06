set -euo pipefail
python - <<'PY'
from pathlib import Path
path=Path('/scratch/result.py')
text=path.read_text()
old = '''    per_donor_cols = [
        "donor_id",
        OUTCOME_COLUMN,
        FEATURE_COLUMN,
        f"{WINNING_VARIATION}__n_reactive_centers",
'''
new = '''    per_donor_cols = [
        "donor_id",
        OUTCOME_COLUMN,
        f"{WINNING_VARIATION}__n_reactive_centers",
'''
text=text.replace(old,new)
path.write_text(text)
PY
python /scratch/result.py | tee /scratch/run_output.txt | sed -n '1,18p'
