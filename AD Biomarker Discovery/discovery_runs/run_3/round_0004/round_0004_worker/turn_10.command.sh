set -euo pipefail
python - <<'PY'
from pathlib import Path
for p in list(Path('/data/autoresearch_runs').rglob('report.md'))[:20]:
    txt=p.read_text(errors='ignore')
    if 'Adjusted Score' in txt or 'adjusted score' in txt:
        print('---',p)
        for line in txt.splitlines():
            if 'Adjusted' in line or 'adjusted' in line:
                print(line)
        break
PY
