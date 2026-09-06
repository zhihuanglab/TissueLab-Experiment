set -euo pipefail
python - <<'PY'
from pathlib import Path
p=Path('/scratch/result.py')
s=p.read_text()
old='''    merged = usable.merge(per_donor, on="donor_id", how="left")\n    residual = merged["outcome"] - merged["predicted"]\n'''
new='''    per_donor_base = per_donor.drop(columns=[spec.feature_column], errors="ignore")\n    merged = usable.merge(per_donor_base, on="donor_id", how="left")\n    residual = merged["outcome"] - merged["predicted"]\n'''
s=s.replace(old,new)
p.write_text(s)
PY
python -m py_compile /scratch/result.py
python /scratch/result.py | tee /scratch/result_stdout.txt
