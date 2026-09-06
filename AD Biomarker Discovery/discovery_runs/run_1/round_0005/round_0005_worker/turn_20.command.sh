python - <<'PY'
from pathlib import Path
path = Path('/scratch/result.py')
text = path.read_text()
old = """    donors = ", ".join(
        f"{row.donor_id} (err={row.abs_error:.3f}, feature={getattr(row, feature_col):.4f})"
        for row in top.itertuples(index=False)
    )
"""
new = """    donors = ", ".join(
        f"{row['donor_id']} (err={row['abs_error']:.3f}, feature={row[feature_col]:.4f})"
        for _, row in top.iterrows()
    )
"""
path.write_text(text.replace(old, new))
PY
python /scratch/result.py
