python - <<'PY'
from pathlib import Path
path = Path('/scratch/result.py')
text = path.read_text()
old = """    winner_table = best_metric[\"loo_table\"].merge(\n        donor_table[[\"donor_id\", \"ca1_pyramidal_fraction\", \"ca1_pyramidal_density\"]],\n        on=\"donor_id\",\n        how=\"left\",\n    )\n"""
new = """    winner_base = best_metric[\"loo_table\"].drop(\n        columns=[col for col in [\"ca1_pyramidal_fraction\", \"ca1_pyramidal_density\"] if col in best_metric[\"loo_table\"].columns],\n        errors=\"ignore\",\n    )\n    winner_table = winner_base.merge(\n        donor_table[[\"donor_id\", \"ca1_pyramidal_fraction\", \"ca1_pyramidal_density\"]],\n        on=\"donor_id\",\n        how=\"left\",\n    )\n"""
text = text.replace(old, new)
path.write_text(text)
PY
python -m py_compile /scratch/result.py
