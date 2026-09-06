set -euo pipefail
python - <<'PY'
from pathlib import Path
path=Path('/scratch/result.py')
text=path.read_text()
text=text.replace(
'''        donor_table[
            [
                "donor_id",
                "slide_name",
                "braak_numeric",
                "cerad_ordinal",
                "cognitive_status",
                "sex",
                "max_age_vis",
                "ec_pyramidal_count",
                "ec_reactive_astrocyte_count",
                f"{feature_col}__proximal_count",
            ]
        ],''',
'''        donor_table[
            [
                c for c in [
                    "donor_id",
                    "slide_name",
                    "braak_numeric",
                    "cerad_ordinal",
                    "sex",
                    "max_age_vis",
                    "ec_pyramidal_count",
                    "ec_reactive_astrocyte_count",
                    f"{feature_col}__proximal_count",
                ] if c in donor_table.columns
            ]
        ],''')
text=text.replace(
'''                f'{row["donor_id"]} (err={row["error"]:+.4f}, '
                f'braak={int(row["braak_numeric"]) if pd.notna(row["braak_numeric"]) else "NA"}, '
                f'cerad={int(row["cerad_ordinal"]) if pd.notna(row["cerad_ordinal"]) else "NA"}, '
                f'{row["cognitive_status"]})' ''',
'''                f'{row["donor_id"]} (err={row["error"]:+.4f}, '
                f'braak={int(row["braak_numeric"]) if pd.notna(row.get("braak_numeric")) else "NA"}, '
                f'cerad={int(row["cerad_ordinal"]) if pd.notna(row.get("cerad_ordinal")) else "NA"}, '
                f'sex={row.get("sex", "NA")})' ''')
text=text.replace(
'''    under_dementia = int((under["cognitive_status"] == "Dementia").sum()) if not under.empty else 0
    over_dementia = int((over["cognitive_status"] == "Dementia").sum()) if not over.empty else 0
    mean_under_prox = float(under[prox_col].mean()) if not under.empty else float("nan")
    mean_over_prox = float(over[prox_col].mean()) if not over.empty else float("nan")
''',
'''    mean_under_prox = float(under[prox_col].mean()) if (not under.empty and prox_col in under.columns) else float("nan")
    mean_over_prox = float(over[prox_col].mean()) if (not over.empty and prox_col in over.columns) else float("nan")
    mean_under_braak = float(under["braak_numeric"].mean()) if (not under.empty and "braak_numeric" in under.columns) else float("nan")
    mean_over_braak = float(over["braak_numeric"].mean()) if (not over.empty and "braak_numeric" in over.columns) else float("nan")
''')
text=text.replace(
'''   Largest underpredictions were {donor_list(under)}. Largest overpredictions were {donor_list(over)}. The underpredicted set had `{under_dementia}`/3 dementia donors versus `{over_dementia}`/3 in the overpredicted set. Mean proximal pyramidal counts were `{mean_under_prox:.1f}` for underpredictions and `{mean_over_prox:.1f}` for overpredictions, suggesting remaining errors are not just sparse-count failures but donor-specific mismatch between EC niche morphology and memory decline.
''',
'''   Largest underpredictions were {donor_list(under)}. Largest overpredictions were {donor_list(over)}. Mean Braak values were `{mean_under_braak:.2f}` for underpredictions and `{mean_over_braak:.2f}` for overpredictions, while mean proximal pyramidal counts were `{mean_under_prox:.1f}` and `{mean_over_prox:.1f}` respectively. That suggests remaining errors are not just sparse-count failures but donor-specific mismatch between EC niche morphology and memory decline.
''')
text=text.replace(
'''        "cognitive_status",''','''        ''')
text=text.replace(
'''            "cognitive_status",''','''            ''')
path.write_text(text)
print('patched')
PY
python /scratch/result.py
echo
echo '--- results.json ---'
sed -n '1,220p' /scratch/results.json
echo
echo '--- report.md ---'
sed -n '1,220p' /scratch/report.md
