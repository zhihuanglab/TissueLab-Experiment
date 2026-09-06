python - <<'PY'
from pathlib import Path
path=Path('/scratch/result.py')
text=path.read_text()
old = '''        donor_table[
            [
                "donor_id",
                feature_col,
                "ca1_pyramidal_n",
                "ca1_reactive_astrocyte_n",
                "ca1_corpora_amylacea_n",
                "reactive_niche_pyramidal_n",
            ]
        ],'''
new = '''        donor_table[
            [
                "donor_id",
                "ca1_pyramidal_n",
                "ca1_reactive_astrocyte_n",
                "ca1_corpora_amylacea_n",
                "reactive_niche_pyramidal_n",
            ]
        ],'''
if old not in text:
    raise SystemExit('pattern not found')
path.write_text(text.replace(old,new))
print('patched merge columns')
PY
