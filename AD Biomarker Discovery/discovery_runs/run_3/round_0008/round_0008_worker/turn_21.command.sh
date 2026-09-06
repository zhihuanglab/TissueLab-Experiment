set -euo pipefail
python - <<'PY'
from pathlib import Path
p=Path('/scratch/result.py')
text=p.read_text()
old="""    predictions_table = best[\"predictions\"].copy()\n    predictions_table = predictions_table.merge(\n        donor_df[[\"donor_id\", best[\"feature_column\"], \"ca1_pyramidal_n\", \"ca1_reactive_astrocyte_n\", \"ca1_lymphocyte_n\"]],\n        on=\"donor_id\",\n        how=\"left\",\n    )\n"""
new="""    predictions_table = best[\"predictions\"].copy()\n    predictions_table = predictions_table.merge(\n        donor_df[[\"donor_id\", \"ca1_pyramidal_n\", \"ca1_reactive_astrocyte_n\", \"ca1_lymphocyte_n\"]],\n        on=\"donor_id\",\n        how=\"left\",\n    )\n"""
if old not in text:
    raise SystemExit('target block not found')
p.write_text(text.replace(old,new))
PY
python /scratch/result.py > /scratch/run_output.txt
sed -n '1,220p' /scratch/run_output.txt
