set -euo pipefail
python - <<'PY'
from pathlib import Path
p = Path('/scratch/result.py')
txt = p.read_text()
old = """    merged = best_loo_table.merge(\n        donor_table.loc[\n            :,\n            [\n                \"donor_id\",\n                \"ca1_astrocyte_count\",\n                \"ca1_reactive_astrocyte_count\",\n                str(best[\"feature_column\"]),\n            ],\n        ],\n        on=\"donor_id\",\n        how=\"left\",\n    )\n"""
new = """    merged = best_loo_table.merge(\n        donor_table.loc[\n            :,\n            [\n                \"donor_id\",\n                \"ca1_astrocyte_count\",\n                \"ca1_reactive_astrocyte_count\",\n            ],\n        ],\n        on=\"donor_id\",\n        how=\"left\",\n    )\n"""
if old not in txt:
    raise SystemExit('merge block not found')
p.write_text(txt.replace(old, new))
PY
python /scratch/result.py > /scratch/run_stdout.txt
sed -n '1,220p' /scratch/report.md
