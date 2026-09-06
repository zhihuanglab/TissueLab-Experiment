set -euo pipefail
python - <<'PY'
from pathlib import Path
p=Path('/scratch/result.py')
text=p.read_text()
# add warnings filter after imports
if 'import warnings' not in text:
    text=text.replace('import re\nimport sys\n', 'import re\nimport sys\nimport warnings\n')
    text=text.replace('import pandas as pd\n\nsys.path.insert(0, "/shared/lib")\n', 'import pandas as pd\n\nwarnings.filterwarnings("ignore", category=UserWarning, module="zarr")\nsys.path.insert(0, "/shared/lib")\n')
# fix direction logic in report generator
old='''    sign = best["partial_r"]\n    if _is_finite(sign) and sign > 0:\n        direction = "Higher"\n        relation = "higher"\n    elif _is_finite(sign) and sign < 0:\n        direction = "Lower values of"\n        relation = "lower"\n    else:\n        direction = "Variation in"\n        relation = "indeterminate"\n'''
if old in text:
    text=text.replace(old,'''    sign = best["partial_r"]\n    if _is_finite(sign) and sign > 0:\n        interpretation_sentence = (\n            f\"Higher values of `{best['feature_column']}` were associated with higher residualized `{OUTCOME_COLUMN}` in this training cohort.\"\n        )\n    elif _is_finite(sign) and sign < 0:\n        interpretation_sentence = (\n            f\"Higher values of `{best['feature_column']}` were associated with lower residualized `{OUTCOME_COLUMN}` in this training cohort.\"\n        )\n    else:\n        interpretation_sentence = (\n            f\"The direction of association between `{best['feature_column']}` and residualized `{OUTCOME_COLUMN}` was indeterminate in this training cohort.\"\n        )\n''')
text=text.replace('{direction} values of `{best[\'feature_column\']}` were associated with {relation} residualized `{OUTCOME_COLUMN}` in this training cohort. ',
                  '{interpretation_sentence} ')
text=text.replace('These errors are {shared_cognition} for cognitive status and {shared_path} for AD neuropath change, so the mistakes appear concentrated in a partly shared clinical/pathology subset rather than being purely random.',
                  'Cohort-level clinical/pathology labels were unavailable in the mounted cohort table, so shared-error interpretation is limited; based on the extracted feature values alone, the misses span both very high and relatively modest CA1 reactive-astrocyte fractions.')
p.write_text(text)

report = Path('/scratch/report.md').read_text()
report = report.replace('These errors are missing for cognitive status and missing for AD neuropath change, so the mistakes appear concentrated in a partly shared clinical/pathology subset rather than being purely random.',
                        'Cohort-level clinical/pathology labels were unavailable in the mounted cohort table, so shared-error interpretation is limited; based on the extracted feature values alone, the misses span both very high and relatively modest CA1 reactive-astrocyte fractions.')
report = report.replace('Lower values of `ca1_reactive_astrocyte_lineage_fraction` were associated with lower residualized `slope_zmem0` in this training cohort.',
                        'Higher values of `ca1_reactive_astrocyte_lineage_fraction` were associated with lower residualized `slope_zmem0` in this training cohort.')
Path('/scratch/report.md').write_text(report)
PY
echo '--- patched report head ---'
sed -n '1,120p' /scratch/report.md
