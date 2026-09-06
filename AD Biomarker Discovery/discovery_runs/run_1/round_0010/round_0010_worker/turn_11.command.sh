set -euo pipefail
python - <<'PY'
import os, glob, json, re
paths=glob.glob('/data/autoresearch_runs/**/*report.md', recursive=True)
for p in paths[:20]:
    txt=open(p,errors='ignore').read()
    if 'Adjusted Score' in txt or 'Selection score' in txt:
        print('PATH',p)
        print('\n'.join(txt.splitlines()[:60]))
        break
PY
