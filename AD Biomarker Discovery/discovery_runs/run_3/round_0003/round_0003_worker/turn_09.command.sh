set -euo pipefail
python - <<'PY'
from pathlib import Path
import tifffile
svs = next(Path('/data').glob('*.svs'))
print('sample svs', svs)
with tifffile.TiffFile(svs) as tf:
    page = tf.pages[0]
    desc = page.description
    print('description head:', desc[:500])
    tags = {tag.name: tag.value for tag in page.tags.values() if tag.name in ['XResolution','YResolution','ResolutionUnit']}
    print('tags', tags)
PY
