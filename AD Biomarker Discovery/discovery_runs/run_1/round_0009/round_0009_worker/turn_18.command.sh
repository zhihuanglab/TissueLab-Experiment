python - <<'PY'
from pathlib import Path
p = Path('/scratch/result.py')
text = p.read_text()
start = text.index('def _rewrite_canonical_source(best_spec: dict[str, object]) -> None:\n')
end = text.index('\n\ndef _write_report', start)
new_func = '''def _rewrite_canonical_source(best_spec: dict[str, object]) -> None:
    source_path = Path(__file__)
    try:
        text = source_path.read_text(encoding="utf-8")
        one_line = json.dumps(best_spec, sort_keys=False)
        text = re.sub(
            r'CANONICAL_VARIATION = .*?# CANONICAL_VARIATION_AUTOFILL',
            f"CANONICAL_VARIATION = {one_line}  # CANONICAL_VARIATION_AUTOFILL",
            text,
        )
        feature_name = str(best_spec["feature_column"])
        text = re.sub(
            r'FEATURE_NAME = ".*?"  # FEATURE_NAME_AUTOFILL',
            f'FEATURE_NAME = "{feature_name}"  # FEATURE_NAME_AUTOFILL',
            text,
        )
        text = re.sub(
            r'FEATURE_COLUMN = ".*?"  # FEATURE_COLUMN_AUTOFILL',
            f'FEATURE_COLUMN = "{feature_name}"  # FEATURE_COLUMN_AUTOFILL',
            text,
        )
        source_path.write_text(text, encoding="utf-8")
    except Exception:
        pass
'''
text = text[:start] + new_func + text[end:]
p.write_text(text)
PY
python -m py_compile /scratch/result.py
