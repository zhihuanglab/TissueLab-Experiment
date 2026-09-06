set -euo pipefail
echo '--- context_bundle.json ---'
python - <<'PY'
import json
from pathlib import Path
p=Path('/scratch/context_bundle.json')
obj=json.loads(p.read_text())
import pprint
pprint.pp(obj)
PY
echo '--- program.md ---'
sed -n '1,220p' /scratch/program.md
echo '--- preseeded result.py ---'
sed -n '1,260p' /scratch/result.py
