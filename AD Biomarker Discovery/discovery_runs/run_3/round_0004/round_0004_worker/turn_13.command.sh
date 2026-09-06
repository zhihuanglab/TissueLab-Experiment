set -euo pipefail
python - <<'PY'
from pathlib import Path
path = Path('/scratch/result.py')
text = path.read_text()
old = """    def donor_summary(frame: pd.DataFrame) -> str:
        if frame.empty:
            return "none"
        return ", ".join(
            f"{row.donor_id} ({row.cognitive_status}, Braak {int(row.braak_numeric)}, {row.sex})"
            for row in frame.itertuples()
        )
"""
new = """    def donor_summary(frame: pd.DataFrame) -> str:
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
if old not in text:
    raise SystemExit('target block not found')
path.write_text(text.replace(old, new))
PY
python /scratch/result.py
python - <<'PY'
from pathlib import Path
print('--- patched canonical lines ---')
for line in Path('/scratch/result.py').read_text().splitlines():
    if line.startswith('CANONICAL_VARIATION = ') or line.startswith('FEATURE_NAME = ') or line.startswith('FEATURE_COLUMN = '):
        print(line)
print('--- wrote files ---')
for p in ['/scratch/results.json','/scratch/report.md','/scratch/donor_feature_table.csv']:
    path=Path(p)
    print(p, path.exists(), path.stat().st_size if path.exists() else None)
PY
