python - <<'PY'
from pathlib import Path
path = Path('/scratch/result.py')
text = path.read_text()
old = "        n_boot=1000,\n        random_state=variation.__hash__() & 0xFFFF,\n    )"
new = "        n_boot=1000,\n        random_state={\"reactive_over_astroglia\": 11, \"reactive_over_all_cells\": 23}[variation],\n    )"
if old not in text:
    raise SystemExit('pattern not found')
path.write_text(text.replace(old,new))
PY
python /scratch/result.py >/scratch/run_stdout.txt
tail -n 20 /scratch/run_stdout.txt
