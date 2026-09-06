set -euo pipefail
python - <<'PY'
from pathlib import Path
p=Path('/scratch/result.py')
s=p.read_text()
old='''    gap = abs(partial_r - loo_r) if np.isfinite(partial_r) and np.isfinite(loo_r) else float("nan")\n'''
new='''    gap = abs(abs(partial_r) - loo_r) if np.isfinite(partial_r) and np.isfinite(loo_r) else float("nan")\n'''
s=s.replace(old,new)
p.write_text(s)
PY
python -m py_compile /scratch/result.py
python /scratch/result.py | tee /scratch/result_stdout.txt
