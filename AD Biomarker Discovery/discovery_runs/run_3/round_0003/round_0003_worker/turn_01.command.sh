set -euo pipefail
echo '--- /scratch/context_bundle.json ---'
python - <<'PY'
import json
from pathlib import Path
p=Path('/scratch/context_bundle.json')
obj=json.loads(p.read_text())
import pprint
pprint.pp(obj)
PY
echo '--- /scratch/program.md ---'
cat /scratch/program.md
echo '--- /scratch listing ---'
ls -lah /scratch
