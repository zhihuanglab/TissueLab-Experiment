set -euo pipefail
python3 - <<'PY'
from pathlib import Path
path = Path('/scratch/result.py')
text = path.read_text()
old = '    gap = float(is_partial_r - loo_predictive_r)\n'
new = '    gap = float(abs(is_partial_r) - abs(loo_predictive_r))\n'
if old not in text:
    raise SystemExit('expected text not found')
path.write_text(text.replace(old, new), encoding='utf-8')
print('patched gap definition')
PY
/usr/local/bin/python3 /scratch/result.py | tee /scratch/run_stdout.txt
grep -n "IS-LOO Gap" -n /scratch/report.md
