set -euo pipefail
python - <<'PY'
from pathlib import Path
path = Path('/scratch/result.py')
text = path.read_text()
text = text.replace('WINNING_VARIATION = "candidate_variant_a"', 'WINNING_VARIATION = "candidate_variant_b"', 1)
text = text.replace('FEATURE_COLUMN = f"{FEATURE_NAME}__{WINNING_VARIATION}"', 'FEATURE_COLUMN = f"{FEATURE_NAME}__{WINNING_VARIATION}"', 1)
text = text.replace(
    '        "The tested depletion contrast produced only weak local single-feature signal, so replacing pyramidal-neuron readout with oligodendrocyte neighbor counts did not recover a cleaner independent axis of injury in this cohort. That suggests the hotspot frame is more tightly aligned with neuron-loss/gliosis than with a separable focal myelin-support deficit."',
    '        "The tested depletion contrast produced only moderate rather than dominant local signal, so replacing pyramidal-neuron readout with oligodendrocyte neighbor counts recovered a real niche effect but not an obviously overwhelming one. That suggests the hotspot frame may capture a myelin-support dimension that is present yet still partly intertwined with the broader CA1 injury program."'
)
text = text.replace(
    'Using the median local oligodendrocyte count is auditable and robust to outliers, but the weak score implies that tightening the hotspot threshold from k=3 to k=4 did not reveal a sharper donor-level effect. The local winner is strong enough that it may contribute information beyond the current CA1 gliosis/neuron-loss panel.',
    'Using the median local oligodendrocyte count is auditable and robust to outliers, and the stronger k=4 hotspot threshold slightly sharpened the donor-level effect versus k=3 while keeping full coverage. The local winner is strong enough that it may contribute information beyond the current CA1 gliosis/neuron-loss panel.'
)
path.write_text(text)
PY
python /scratch/result.py > /scratch/run_output.txt
sed -n '1,22p' /scratch/run_output.txt
grep -n 'WINNING_VARIATION' /scratch/result.py | head
sed -n '1,200p' /scratch/report.md | sed -n '1,80p'
