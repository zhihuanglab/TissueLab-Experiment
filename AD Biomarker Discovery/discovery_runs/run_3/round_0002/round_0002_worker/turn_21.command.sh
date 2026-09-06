python - <<'PY'
from pathlib import Path
path=Path('/scratch/result.py')
text=path.read_text()
text=text.replace(
'        merged.loc[:, ["donor_id", "cognitive_status", "braak_numeric", "cerad_ordinal", "sex", f"{best[\'feature_column\']}_denom"]],\n',
'        merged.loc[:, ["donor_id", "braak_numeric", "cerad_ordinal", "sex", f"{best[\'feature_column\']}_denom"]],\n'
)
text=text.replace(
'    top_records = top.loc[:, ["donor_id", "outcome", "predicted", "abs_error", "cognitive_status", "braak_numeric", "cerad_ordinal", "sex", f"{best[\'feature_column\']}_denom"]].to_dict(orient="records")\n    common_status = top["cognitive_status"].mode(dropna=True)\n    shared_status = str(common_status.iloc[0]) if not common_status.empty else "mixed cognitive-status"\n    median_braak = float(top["braak_numeric"].median()) if len(top) else float("nan")\n    median_denom = float(top[f"{best[\'feature_column\']}_denom"].median()) if len(top) else float("nan")\n    summary = (\n        f"The largest LOO errors cluster in donors with {shared_status} status, "\n        f"median Braak {median_braak:.1f}, and a median peripyramidal denominator of {median_denom:.0f}; "\n        f"these appear to be donors where CA1 reactive fraction alone does not fully track outcome severity."\n    )\n',
'    top_records = top.loc[:, ["donor_id", "outcome", "predicted", "abs_error", "braak_numeric", "cerad_ordinal", "sex", f"{best[\'feature_column\']}_denom"]].to_dict(orient="records")\n    sex_mode = top["sex"].mode(dropna=True)\n    shared_sex = str(sex_mode.iloc[0]) if not sex_mode.empty else "mixed-sex"\n    median_braak = float(top["braak_numeric"].median()) if len(top) else float("nan")\n    median_denom = float(top[f"{best[\'feature_column\']}_denom"].median()) if len(top) else float("nan")\n    summary = (\n        f"The largest LOO errors are concentrated in {shared_sex} donors with "\n        f"median Braak {median_braak:.1f} and a median peripyramidal denominator of {median_denom:.0f}; "\n        f"these appear to be donors where CA1 reactive fraction alone does not fully track outcome severity."\n    )\n'
)
text=text.replace(
'            f"|error| {row[\'abs_error\']:.3f}, status {row[\'cognitive_status\']}, Braak {row[\'braak_numeric\']}, "\n            f"CERAD {row[\'cerad_ordinal\']}, denominator {int(row[denom_key])}."\n',
'            f"|error| {row[\'abs_error\']:.3f}, sex {row[\'sex\']}, Braak {row[\'braak_numeric\']}, "\n            f"CERAD {row[\'cerad_ordinal\']}, denominator {int(row[denom_key])}."\n'
)
path.write_text(text)
PY
python -m py_compile /scratch/result.py
python /scratch/result.py > /scratch/run_stdout.txt
cat /scratch/run_stdout.txt
echo '--- report.md head ---'
sed -n '1,220p' /scratch/report.md
