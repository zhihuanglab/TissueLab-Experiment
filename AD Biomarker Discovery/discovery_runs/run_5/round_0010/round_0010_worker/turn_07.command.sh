set -euo pipefail
python - <<'PY'
import tifffile
from pathlib import Path
p=next(Path('/data').glob('*.svs'))
with tifffile.TiffFile(str(p)) as tf:
    desc=tf.pages[0].description or ''
    print(p.name)
    print(desc[:500])
PY
