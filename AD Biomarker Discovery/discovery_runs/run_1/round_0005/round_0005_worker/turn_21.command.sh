python - <<'PY'
from pathlib import Path
path=Path('/scratch/result.py')
text=path.read_text()
text=text.replace(
"""        donor_table[
            [
                "donor_id",
                feature_col,
                "ca1_pyramidal_count",
                "ca1_reactive_astrocyte_count",
            ]
        ],
""",
"""        donor_table[
            [
                "donor_id",
                "ca1_pyramidal_count",
                "ca1_reactive_astrocyte_count",
            ]
        ],
""")
path.write_text(text)
PY
python /scratch/result.py
