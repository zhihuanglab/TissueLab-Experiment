set -euo pipefail
python - <<'PY'
from pathlib import Path
txt=Path('/shared/templates/worker_embedding_mechanistic_template.py').read_text().splitlines()
for start in [70,80,160,180,194]:
    print(f'--- lines {start}-{start+35} ---')
    for i in range(start-1,min(start+35,len(txt))):
        print(f'{i+1}: {txt[i]}')
PY
