python - <<'PY'
from pathlib import Path
p = Path('/scratch/result.py')
text = p.read_text()
old = """    for row in cohort.itertuples(index=False):
        slide_path = data_root / str(row.slide_name)
        features = _extract_slide_family_features(slide_path)
        rows.append(
            {
                "donor_id": str(row.donor_id),
                "slide_name": str(row.slide_name),
                "cognitive_status": str(row.cognitive_status),
                "sex": str(row.sex),
                "sex_binary": float(row.sex_binary) if pd.notna(row.sex_binary) else float("nan"),
                OUTCOME_COLUMN: float(getattr(row, OUTCOME_COLUMN)),
                "max_age_vis": float(row.max_age_vis),
                "braak_numeric": float(row.braak_numeric),
                "cerad_ordinal": float(row.cerad_ordinal),
                **features,
            }
        )
"""
new = """    for _, row in cohort.iterrows():
        slide_path = data_root / str(row["slide_name"])
        features = _extract_slide_family_features(slide_path)
        rows.append(
            {
                "donor_id": str(row["donor_id"]),
                "slide_name": str(row["slide_name"]),
                "cognitive_status": str(row["cognitive_status"]),
                "sex": str(row["sex"]),
                "sex_binary": float(row["sex_binary"]) if pd.notna(row["sex_binary"]) else float("nan"),
                OUTCOME_COLUMN: float(row[OUTCOME_COLUMN]),
                "max_age_vis": float(row["max_age_vis"]),
                "braak_numeric": float(row["braak_numeric"]),
                "cerad_ordinal": float(row["cerad_ordinal"]),
                **features,
            }
        )
"""
if old not in text:
    raise SystemExit('target block not found')
p.write_text(text.replace(old,new))
PY
python /scratch/result.py
