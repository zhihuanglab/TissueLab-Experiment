python - <<'PY'
from pathlib import Path
p=Path('/scratch/result.py')
text=p.read_text()
text=text.replace('    frame = df.loc[:, ["donor_id", "slide_name", "cognitive_status", outcome_col, feature_col, *confounds]].dropna().reset_index(drop=True)\n',
                  '    keep_cols = ["donor_id", "slide_name", outcome_col, feature_col, *confounds]\n    if "cognitive_status" in df.columns:\n        keep_cols.insert(2, "cognitive_status")\n    frame = df.loc[:, keep_cols].dropna().reset_index(drop=True)\n')
text=text.replace('                "cognitive_status": str(test.iloc[0]["cognitive_status"]),\n', '')
text=text.replace('                "cognitive_status": str(row["cognitive_status"]),\n', '')
text=text.replace('        counts = table.loc[table["donor_id"].isin(top_err["donor_id"]), ["donor_id", "cognitive_status", "n_ca1_pyramidal_proximal_ra_80px"]].copy()\n',
                  '        counts = table.loc[table["donor_id"].isin(top_err["donor_id"]), ["donor_id", "n_ca1_pyramidal_proximal_ra_80px"]].copy()\n')
text=text.replace('        "cognitive_status",\n', '')
text=text.replace('        "cognitive_status",\n        "sex",\n', '        "sex",\n')
p.write_text(text)
PY
python /scratch/result.py
