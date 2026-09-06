python - <<'PY'
from pathlib import Path
path=Path('/scratch/result.py')
text=path.read_text()
text=text.replace(
    "    for row in top_errors:\n        top_error_lines.append(\n            f\"- {row['donor_id']}: outcome {row['outcome']:.3f}, predicted {row['predicted']:.3f}, \"\n            f\"|error| {row['abs_error']:.3f}, status {row['cognitive_status']}, Braak {row['braak_numeric']}, \"\n            f\"CERAD {row['cerad_ordinal']}, denominator {int(row[f\\\"{best['feature_column']}_denom\\\"])}.\"\n        )\n",
    "    denom_key = f\"{best['feature_column']}_denom\"\n    for row in top_errors:\n        top_error_lines.append(\n            f\"- {row['donor_id']}: outcome {row['outcome']:.3f}, predicted {row['predicted']:.3f}, \"\n            f\"|error| {row['abs_error']:.3f}, status {row['cognitive_status']}, Braak {row['braak_numeric']}, \"\n            f\"CERAD {row['cerad_ordinal']}, denominator {int(row[denom_key])}.\"\n        )\n"
)
path.write_text(text)
PY
python -m py_compile /scratch/result.py
