python - <<'PY'
from pathlib import Path
path = Path('/scratch/result.py')
text = path.read_text()

old = """        merge_cols = [\"donor_id\", \"cognitive_status\", \"braak_label\", \"cerad_label\", \"overall_ad_neuropath_change\", \"sex\"]\n        merged = winner_table.merge(donor_table[merge_cols], on=\"donor_id\", how=\"left\")\n"""
new = """        merge_cols = [\n            col\n            for col in [\"donor_id\", \"cognitive_status\", \"braak_label\", \"cerad_label\", \"overall_ad_neuropath_change\", \"sex\"]\n            if col in donor_table.columns\n        ]\n        merged = winner_table.merge(donor_table[merge_cols], on=\"donor_id\", how=\"left\")\n"""
text = text.replace(old, new)

old2 = """        for col, label in [\n            (\"cognitive_status\", \"cognitive status\"),\n            (\"overall_ad_neuropath_change\", \"AD neuropath change\"),\n            (\"braak_label\", \"Braak\"),\n            (\"cerad_label\", \"CERAD\"),\n            (\"sex\", \"sex\"),\n        ]:\n"""
new2 = """        for col, label in [\n            (\"cognitive_status\", \"cognitive status\"),\n            (\"overall_ad_neuropath_change\", \"AD neuropath change\"),\n            (\"braak_label\", \"Braak\"),\n            (\"cerad_label\", \"CERAD\"),\n            (\"sex\", \"sex\"),\n        ]:\n"""
# no-op, keeps block for clarity
text = text.replace(old2, new2)

path.write_text(text)
PY
python -m py_compile /scratch/result.py
