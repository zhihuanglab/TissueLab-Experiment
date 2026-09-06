python - <<'PY'
from pathlib import Path
path=Path('/scratch/result.py')
text=path.read_text()
text=text.replace(
'        reactive = int(((astro_lineage["cell_type"] == "Reactive Astrocyte").to_numpy()) & near_mask).sum() if denom else 0\n',
'        reactive = int((((astro_lineage["cell_type"] == "Reactive Astrocyte").to_numpy()) & near_mask).sum()) if denom else 0\n'
)
# suppress repetitive zarr warnings
if 'import warnings' not in text:
    text=text.replace('import json\n', 'import json\nimport warnings\n')
    text=text.replace('from shared_analysis import build_cell_table, load_training_cohort\n', 'from shared_analysis import build_cell_table, load_training_cohort\n\nwarnings.filterwarnings("ignore", category=UserWarning, module="zarr")\n')
path.write_text(text)
PY
python -m py_compile /scratch/result.py
python /scratch/result.py > /scratch/run_stdout.txt
cat /scratch/run_stdout.txt
echo '--- results.json head ---'
python - <<'PY'
import json
from pathlib import Path
p=Path('/scratch/results.json')
print('exists', p.exists(), 'size', p.stat().st_size if p.exists() else None)
if p.exists():
    obj=json.loads(p.read_text())
    print({k: obj.get(k) for k in ['best_variation','feature_column','partial_r','selection_score','loo_predictive_r','is_loo_gap','adjusted_score']})
    print('ranked count', len(obj.get('ranked_variations', [])))
PY
echo '--- report.md head ---'
sed -n '1,220p' /scratch/report.md
