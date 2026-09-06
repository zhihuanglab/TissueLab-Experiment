python - <<'PY'
from pathlib import Path
path = Path('/scratch/result.py')
text = path.read_text()
old = """        print(
            f"  {row['donor_id']}  {row['outcome']:.4f}  "
            f"{row['predicted']:.4f if math.isfinite(float(row['predicted'])) else 'nan'}  "
            f"{row[best_eval.feature_column]:.6f}"
        )
"""
new = """        pred = row['predicted']
        pred_text = f"{float(pred):.4f}" if pred is not None and math.isfinite(float(pred)) else "nan"
        print(
            f"  {row['donor_id']}  {row['outcome']:.4f}  "
            f"{pred_text}  "
            f"{row[best_eval.feature_column]:.6f}"
        )
"""
if old not in text:
    raise SystemExit('target block not found')
path.write_text(text.replace(old,new))
PY
python -m py_compile /scratch/result.py
