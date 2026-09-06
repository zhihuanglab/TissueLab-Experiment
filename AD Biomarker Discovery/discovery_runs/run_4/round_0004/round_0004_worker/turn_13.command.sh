set -euo pipefail
python - <<'PY'
from pathlib import Path
p = Path('/scratch/result.py')
txt = p.read_text()
old = """        for row in err_rows.itertuples(index=False):\n            astro_n = int(getattr(row, \"ca1_astrocyte_count\", float(\"nan\")) or 0)\n            react_n = int(getattr(row, \"ca1_reactive_astrocyte_count\", float(\"nan\")) or 0)\n            if min(astro_n, react_n) < 2 * MIN_CELLS_PER_GROUP:\n                near_threshold += 1\n            donor_bits.append(\n                f\"{row.donor_id} (obs={getattr(row, OUTCOME_COLUMN):.3f}, pred={row.predicted:.3f}, \"\n                f\"feature={getattr(row, str(best['feature_column'])):.3f}, astro={astro_n}, reactive={react_n})\"\n            )\n"""
new = """        for row in err_rows.itertuples(index=False):\n            row_dict = row._asdict()\n            astro_n = int(float(row_dict.get(\"ca1_astrocyte_count\", 0.0) or 0.0))\n            react_n = int(float(row_dict.get(\"ca1_reactive_astrocyte_count\", 0.0) or 0.0))\n            if min(astro_n, react_n) < 2 * MIN_CELLS_PER_GROUP:\n                near_threshold += 1\n            donor_bits.append(\n                f\"{row_dict.get('donor_id')} (obs={float(row_dict.get(OUTCOME_COLUMN, float('nan'))):.3f}, \"\n                f\"pred={float(row_dict.get('predicted', float('nan'))):.3f}, \"\n                f\"feature={float(row_dict.get(str(best['feature_column']), float('nan'))):.3f}, \"\n                f\"astro={astro_n}, reactive={react_n})\"\n            )\n"""
if old not in txt:
    raise SystemExit('target block not found')
p.write_text(txt.replace(old, new))
PY
python /scratch/result.py
