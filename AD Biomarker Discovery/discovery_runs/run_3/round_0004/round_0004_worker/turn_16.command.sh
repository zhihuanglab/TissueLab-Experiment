set -euo pipefail
python - <<'PY'
from pathlib import Path
path = Path('/scratch/result.py')
text = path.read_text()
old = """    def donor_summary(frame: pd.DataFrame) -> str:
        if frame.empty:
            return "none"
        parts = []
        for _, row in frame.iterrows():
            cognition = row.get("cognitive_status", "NA")
            braak = row.get("braak_numeric", np.nan)
            sex = row.get("sex", "NA")
            braak_text = "NA" if pd.isna(braak) else str(int(braak))
            parts.append(f"{row['donor_id']} ({cognition}, Braak {braak_text}, {sex})")
        return ", ".join(parts)
"""
new = """    def donor_summary(frame: pd.DataFrame) -> str:
        if frame.empty:
            return "none"
        parts = []
        for _, row in frame.iterrows():
            descriptors = []
            if "cognitive_status" in frame.columns and pd.notna(row.get("cognitive_status", np.nan)):
                descriptors.append(str(row["cognitive_status"]))
            if "braak_numeric" in frame.columns and pd.notna(row.get("braak_numeric", np.nan)):
                descriptors.append(f"Braak {int(row['braak_numeric'])}")
            if "sex" in frame.columns and pd.notna(row.get("sex", np.nan)):
                descriptors.append(str(row["sex"]))
            desc = ", ".join(descriptors)
            parts.append(f"{row['donor_id']} ({desc})" if desc else str(row["donor_id"]))
        return ", ".join(parts)
"""
if old not in text:
    raise SystemExit('target block not found for donor_summary')
path.write_text(text.replace(old, new))
PY
python /scratch/result.py >/scratch/run_stdout.txt
sed -n '1,30p' /scratch/run_stdout.txt
echo '--- report excerpt ---'
sed -n '1,120p' /scratch/report.md
