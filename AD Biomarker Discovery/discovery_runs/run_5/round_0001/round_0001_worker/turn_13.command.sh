set -euo pipefail
python - <<'PY'
import sys
sys.path.insert(0,'/shared/lib')
from shared_analysis.sea_ad_lfb import load_training_cohort
from pathlib import Path
cohort=load_training_cohort(Path('/data'))
print(cohort.columns.tolist())
PY
python - <<'PY'
from pathlib import Path
p=Path('/scratch/result.py')
s=p.read_text()
s=s.replace('''            "sex": row["sex"],\n            "cognitive_status": row["cognitive_status"],\n''','''            "sex": row["sex"],\n''')
p.write_text(s)
PY
python -m py_compile /scratch/result.py
python /scratch/result.py | tee /scratch/result_stdout.txt
