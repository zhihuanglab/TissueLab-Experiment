python - <<'PY'
from pathlib import Path
path = Path('/scratch/result.py')
text = path.read_text(encoding='utf-8')
text = text.replace(
'    load_class_lookup,\n',
'    load_class_ids,\n    load_class_lookup,\n',
)
text = text.replace(
'''    centroids = load_centroids(slide_path)
    class_ids = load_class_lookup(slide_path)
    name_to_id = {name: idx for idx, name in class_ids.items()}
    pyramidal_id = name_to_id.get(PYRAMIDAL_CELL_TYPE, None)
    reactive_id = name_to_id.get(REACTIVE_CELL_TYPE, None)
    nuc_ids = np.asarray(
        pd.Series(load_class_lookup(slide_path)).index
    )  # placeholder to satisfy linter; overwritten below
    del nuc_ids
    from shared_analysis import load_class_ids  # local import to keep top-level import short

    nuclei_class_ids = load_class_ids(slide_path)
''',
'''    centroids = load_centroids(slide_path)
    class_lookup = load_class_lookup(slide_path)
    name_to_id = {name: idx for idx, name in class_lookup.items()}
    pyramidal_id = name_to_id.get(PYRAMIDAL_CELL_TYPE, None)
    reactive_id = name_to_id.get(REACTIVE_CELL_TYPE, None)
    nuclei_class_ids = load_class_ids(slide_path)
''')
text = text.replace(
'            random_state=variation_name.__hash__() % (2**32),\n',
'            random_state=(101 if variation_name == "candidate_variant_a" else 202),\n',
)
path.write_text(text, encoding='utf-8')
PY
python /scratch/result.py > /scratch/stdout.txt
cat /scratch/stdout.txt
