set -euo pipefail
python - <<'PY'
from pathlib import Path
path = Path('/scratch/result.py')
text = path.read_text()
old = """    is_loo_gap = None
    if partial_r is not None and _safe_float(loo_r) is not None:
        is_loo_gap = partial_r - float(loo_r)
"""
new = """    is_loo_gap = None
    if partial_r is not None and _safe_float(loo_r) is not None:
        is_loo_gap = abs(partial_r) - float(loo_r)
"""
if old not in text:
    raise SystemExit('target block not found')
path.write_text(text.replace(old, new))
PY

python /scratch/result.py | tee /scratch/result_stdout.txt

python - <<'PY'
import json, re
from pathlib import Path
result_path = Path('/scratch/result.py')
results = json.loads(Path('/scratch/results.json').read_text())
best = results['best_variation']
text = result_path.read_text()
text = re.sub(r'CANONICAL_REPLAY_VARIATION = ".*?"', f'CANONICAL_REPLAY_VARIATION = "{best}"', text, count=1)
result_path.write_text(text)
print(f"Updated canonical replay variation to {best}")
PY
