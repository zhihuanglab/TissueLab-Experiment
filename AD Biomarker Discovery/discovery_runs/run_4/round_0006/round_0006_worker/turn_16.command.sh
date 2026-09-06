python - <<'PY'
from pathlib import Path
p=Path('/scratch/result.py')
txt=p.read_text()
old='''        error_df = loo_frame.merge(
            df[["donor_id", "braak_numeric", "cerad_ordinal", "cognitive_status", "sex"]].drop_duplicates(),
            on="donor_id",
            how="left",
        ).sort_values("abs_error", ascending=False)
'''
new='''        meta_cols = ["donor_id", "braak_numeric", "cerad_ordinal", "sex"]
        if "cognitive_status" in df.columns:
            meta_cols.append("cognitive_status")
        error_df = loo_frame.merge(
            df[meta_cols].drop_duplicates(),
            on="donor_id",
            how="left",
        ).sort_values("abs_error", ascending=False)
'''
txt=txt.replace(old,new)
txt=txt.replace("status {row.get('cognitive_status')}.", "status {row.get('cognitive_status', 'NA')}.")
p.write_text(txt)
PY
python /scratch/result.py
