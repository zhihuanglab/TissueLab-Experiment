set -e
python - <<'PY'
import pandas as pd
cohort = pd.read_csv('/data/training_cohort.csv')
print(list(cohort.columns))
PY
# patch result.py to tolerate missing cognitive_status
python - <<'PY'
from pathlib import Path
p=Path('/scratch/result.py')
text=p.read_text()
old='''def summarize_error_pattern(best_per_donor: list[dict], cohort: pd.DataFrame) -> str:
    if not best_per_donor:
        return "No analyzable donors."
    err = pd.DataFrame(best_per_donor).sort_values("abs_error", ascending=False).head(3)
    merged = err.merge(
        cohort[["donor_id", "cognitive_status", "braak_numeric", "cerad_ordinal", "sex"]],
        on="donor_id",
        how="left",
    )
    donors = ", ".join(merged["donor_id"].tolist())
    dementia_n = int((merged["cognitive_status"] == "Dementia").sum())
    high_braak_n = int((merged["braak_numeric"] >= 5).sum())
    return (
        f"Largest absolute LOO errors were {donors}; "
        f"{dementia_n}/3 are dementia donors and {high_braak_n}/3 have Braak >=5."
    )
'''
new='''def summarize_error_pattern(best_per_donor: list[dict], cohort: pd.DataFrame) -> str:
    if not best_per_donor:
        return "No analyzable donors."
    err = pd.DataFrame(best_per_donor).sort_values("abs_error", ascending=False).head(3)
    extra_cols = [c for c in ["donor_id", "cognitive_status", "braak_numeric", "cerad_ordinal", "sex"] if c in cohort.columns]
    merged = err.merge(cohort[extra_cols], on="donor_id", how="left")
    donors = ", ".join(merged["donor_id"].tolist())
    if "cognitive_status" in merged.columns:
        dementia_n = int((merged["cognitive_status"] == "Dementia").sum())
        dementia_part = f"{dementia_n}/3 are dementia donors and "
    else:
        dementia_part = ""
    if "braak_numeric" in merged.columns:
        high_braak_n = int((merged["braak_numeric"] >= 5).sum())
        braak_part = f"{high_braak_n}/3 have Braak >=5."
    else:
        braak_part = "their shared pathology labels are not all available in the cohort table."
    return f"Largest absolute LOO errors were {donors}; {dementia_part}{braak_part}"
'''
if old not in text:
    raise SystemExit('Expected block not found')
p.write_text(text.replace(old,new))
PY
python -u /scratch/result.py
echo
echo '--- results.json ---'
python - <<'PY'
import json
obj=json.load(open('/scratch/results.json'))
print(json.dumps({k: obj[k] for k in ['best_variation','feature_column','partial_r','selection_score','loo_predictive_r','is_loo_gap','penalty','adjusted_score']}, indent=2))
PY
echo
echo '--- report.md ---'
sed -n '1,220p' /scratch/report.md
