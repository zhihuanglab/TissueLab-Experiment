set -euo pipefail
python - <<'PY'
import sys
sys.path.append('/shared/lib')
from shared_analysis import load_training_cohort
from pathlib import Path
cohort=load_training_cohort(Path('/data'))
print(cohort.columns.tolist())
PY
python - <<'PY'
from pathlib import Path
p=Path('/scratch/result.py')
text=p.read_text()
old='''def _error_pattern_text(loo_table: pd.DataFrame, donor_table: pd.DataFrame) -> str:
    if loo_table.empty:
        return "LOO evaluation was unavailable because too few analyzable donors remained after confound filtering."
    merged = loo_table.merge(
        donor_table.loc[:, ["donor_id", "cognitive_status", "braak_numeric", "cerad_ordinal", "sex"]],
        on="donor_id",
        how="left",
    )
'''
new='''def _error_pattern_text(loo_table: pd.DataFrame, donor_table: pd.DataFrame) -> str:
    if loo_table.empty:
        return "LOO evaluation was unavailable because too few analyzable donors remained after confound filtering."
    meta_cols = ["donor_id"]
    for col in ["cognitive_status", "braak_numeric", "cerad_ordinal", "sex"]:
        if col in donor_table.columns:
            meta_cols.append(col)
    merged = loo_table.merge(
        donor_table.loc[:, meta_cols],
        on="donor_id",
        how="left",
    )
'''
if old not in text:
    raise SystemExit('target block not found')
text=text.replace(old,new)
old='''        donor_bits.append(
            f"{row.donor_id} (abs err={row.abs_error:.3f}, status={row.cognitive_status}, Braak={row.braak_numeric}, CERAD={row.cerad_ordinal})"
        )
    if top.empty:
        return "No donor-level LOO predictions were available."
    status_counts = top["cognitive_status"].fillna("unknown").value_counts().to_dict()
'''
new='''        status = getattr(row, "cognitive_status", "unknown")
        braak = getattr(row, "braak_numeric", "unknown")
        cerad = getattr(row, "cerad_ordinal", "unknown")
        donor_bits.append(
            f"{row.donor_id} (abs err={row.abs_error:.3f}, status={status}, Braak={braak}, CERAD={cerad})"
        )
    if top.empty:
        return "No donor-level LOO predictions were available."
    if "cognitive_status" in top.columns:
        status_counts = top["cognitive_status"].fillna("unknown").value_counts().to_dict()
    else:
        status_counts = {}
'''
if old not in text:
    raise SystemExit('second target block not found')
text=text.replace(old,new)
old='''    if len(status_counts) == 1:
        shared_status = next(iter(status_counts))
        shared_text = f"The largest errors were concentrated in donors labeled {shared_status}."
    else:
        shared_text = "The largest errors were mixed across cognitive-status groups rather than confined to one label."
'''
new='''    if len(status_counts) == 1:
        shared_status = next(iter(status_counts))
        shared_text = f"The largest errors were concentrated in donors labeled {shared_status}."
    elif len(status_counts) > 1:
        shared_text = "The largest errors were mixed across cognitive-status groups rather than confined to one label."
    else:
        shared_text = "The largest errors did not map cleanly to a single available metadata label."
'''
if old not in text:
    raise SystemExit('third target block not found')
text=text.replace(old,new)
p.write_text(text)
print('patched result.py')
PY
python /scratch/result.py
echo '--- results.json head ---'
python - <<'PY'
import json
with open('/scratch/results.json') as f:
    obj=json.load(f)
print(json.dumps({
    'best_variation': obj.get('best_variation'),
    'feature_column': obj.get('feature_column'),
    'selection_score': obj.get('selection_score'),
    'ranked_variations': obj.get('ranked_variations'),
}, indent=2)[:4000])
PY
echo '--- report.md ---'
sed -n '1,220p' /scratch/report.md
