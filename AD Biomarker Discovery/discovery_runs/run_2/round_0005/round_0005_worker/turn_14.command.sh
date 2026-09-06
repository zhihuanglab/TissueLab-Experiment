python - <<'PY'
from pathlib import Path
p = Path('/scratch/result.py')
text = p.read_text()
text = text.replace(
"""        merged = top_err.merge(
            full_table[
                [
                    "donor_id",
                    "cognitive_status",
                    "braak_numeric",
                    "cerad_ordinal",
                    "ca1_reactive_astrocyte_count",
                    "ca1_astrocyte_count",
                    reactive_col,
                    astro_col,
                    best_feature_col,
                ]
            ],
            on="donor_id",
            how="left",
            suffixes=("", "_feat"),
        )""",
"""        merged = top_err.merge(
            full_table[
                [
                    "donor_id",
                    "braak_numeric",
                    "cerad_ordinal",
                    "ca1_reactive_astrocyte_count",
                    "ca1_astrocyte_count",
                    reactive_col,
                    astro_col,
                    best_feature_col,
                ]
            ],
            on="donor_id",
            how="left",
            suffixes=("", "_feat"),
        )"""
)
text = text.replace(
"""                f"{row['donor_id']} (error={format_float(row['abs_error'])}, braak={row['braak_numeric']}, CERAD={row['cerad_ordinal']}, status={row['cognitive_status']}, reactive_frac={format_float(row[reactive_col])}, astro_frac={format_float(row[astro_col])})"
""",
"""                f"{row['donor_id']} (error={format_float(row['abs_error'])}, braak={row['braak_numeric']}, CERAD={row['cerad_ordinal']}, reactive_frac={format_float(row[reactive_col])}, astro_frac={format_float(row[astro_col])})"
"""
)
text = text.replace(
"""        extra_columns=[
            "cognitive_status",
            "ca1_cell_count",""",
"""        extra_columns=[
            "ca1_cell_count","""
)
p.write_text(text)
print('patched')
PY
python -m py_compile /scratch/result.py
