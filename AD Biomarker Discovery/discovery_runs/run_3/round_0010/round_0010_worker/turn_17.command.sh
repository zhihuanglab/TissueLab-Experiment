set -euo pipefail
python - <<'PY'
from pathlib import Path
p=Path('/scratch/result.py')
txt=p.read_text()
old = """    meta_cols = [
        \"donor_id\",
        \"cognitive_status\",
        \"overall_ad_neuropath_change\",
        \"braak_numeric\",
        \"cerad_ordinal\",
        \"ca1_severe_reactive_cuffed_pyramidal_n_r50um_ra2\",
        \"ca1_conditional_immune_cuff_positive_n_r50um_ra2_ly2\",
        \"ca1_conditional_immune_cuff_positive_n_r50um_ra2_ly3\",
    ]
"""
new = """    meta_cols = [
        \"donor_id\",
        \"cognitive_status\",
        \"overall_ad_neuropath_change\",
        \"braak_numeric\",
        \"cerad_ordinal\",
    ]
"""
txt=txt.replace(old,new)
p.write_text(txt)
print('patched2')
PY

python /scratch/result.py > /scratch/run_stdout.txt
printf '%s\n' '--- stdout ---'
sed -n '1,220p' /scratch/run_stdout.txt
printf '\n--- results.json head ---\n'
sed -n '1,180p' /scratch/results.json
printf '\n--- report.md ---\n'
sed -n '1,220p' /scratch/report.md
