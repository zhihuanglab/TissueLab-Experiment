python - <<'PY'
from pathlib import Path
path=Path('/scratch/result.py')
text=path.read_text()
text=text.replace('BEST_VARIATION = "reactive_purity_r40um"\nFEATURE_NAME = "ca1_reactive_local_purity"\nFEATURE_COLUMN = "ca1_reactive_local_purity_r40um"\n',
                  'BEST_VARIATION = "reactive_purity_r60um"\nFEATURE_NAME = "ca1_reactive_local_purity"\nFEATURE_COLUMN = "ca1_reactive_local_purity_r60um"\n')
path.write_text(text)
PY
python /scratch/result.py > /scratch/result_stdout.txt
sed -n '1,60p' /scratch/result_stdout.txt
