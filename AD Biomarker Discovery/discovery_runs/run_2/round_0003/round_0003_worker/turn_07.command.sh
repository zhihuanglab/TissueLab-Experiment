set -euo pipefail
python - <<'PY'
from pathlib import Path
p=Path('/shared/templates/worker_embedding_mechanistic_template.py')
text=p.read_text()
for start,end in [(70,115),(175,220)]:
    print(f'--- lines {start}-{end} ---')
    for i,line in enumerate(text.splitlines()[start-1:end], start):
        print(f'{i:03d}: {line}')
PY
