python - <<'PY'
from pathlib import Path
path = Path('/scratch/result.py')
text = path.read_text()
old = '''## Next
Run one more local sweep in the same CA1 pyramidal–reactive astrocyte niche, keeping the 25th-percentile tail definition fixed while varying the proximity radius around 60–100 px to test whether the signal is truly immediate-neighbor atrophy or a broader peri-reactive-astrocyte field effect.
'''
new = '''## Next
Run one more local sweep in the same CA1 pyramidal–reactive astrocyte niche, keeping the winning tail definition fixed while varying the proximity radius around 60–100 px to test whether the signal is truly immediate-neighbor atrophy or a broader peri-reactive-astrocyte field effect.
'''
if old not in text:
    raise SystemExit('target next block not found')
path.write_text(text.replace(old, new))
PY
python -m py_compile /scratch/result.py
python - <<'PY'
from pathlib import Path
report = Path('/scratch/report.md').read_text()
report = report.replace(
"Run one more local sweep in the same CA1 pyramidal–reactive astrocyte niche, keeping the 25th-percentile tail definition fixed while varying the proximity radius around 60–100 px to test whether the signal is truly immediate-neighbor atrophy or a broader peri-reactive-astrocyte field effect.",
"Run one more local sweep in the same CA1 pyramidal–reactive astrocyte niche, keeping the winning tail definition fixed while varying the proximity radius around 60–100 px to test whether the signal is truly immediate-neighbor atrophy or a broader peri-reactive-astrocyte field effect."
)
Path('/scratch/report.md').write_text(report)
print(Path('/scratch/report.md').read_text())
PY
