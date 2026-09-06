set -euo pipefail
python - <<'PY'
from pathlib import Path
path = Path('/scratch/result.py')
text = path.read_text()
text = text.replace(
"""    loo_df = loo_df.merge(
        table.loc[
            :,
            [
                "donor_id",
                feature_col,
                spec["reactive_summary_column"],
                spec["astrocyte_summary_column"],
                "ca1_reactive_astrocyte_count",
                "ca1_astrocyte_count",
                "cognitive_status",
                "overall_ad_neuropath_change",
            ],
        ],
        on="donor_id",
        how="left",
    )
""",
"""    merge_columns = [
        "donor_id",
        feature_col,
        spec["reactive_summary_column"],
        spec["astrocyte_summary_column"],
        "ca1_reactive_astrocyte_count",
        "ca1_astrocyte_count",
    ]
    for optional_col in ["sex", "cognitive_status", "overall_ad_neuropath_change"]:
        if optional_col in table.columns:
            merge_columns.append(optional_col)
    loo_df = loo_df.merge(
        table.loc[:, merge_columns],
        on="donor_id",
        how="left",
    )
""")
text = text.replace(
"""                "ca1_reactive_astrocyte_plus_astrocyte_count",
                "cognitive_status",
                "overall_ad_neuropath_change",
""",
"""                "ca1_reactive_astrocyte_plus_astrocyte_count",
                "sex",
""")
path.write_text(text)
PY

python /scratch/result.py | tee /scratch/result_stdout.txt

python - <<'PY'
import json
from pathlib import Path
result_path = Path('/scratch/result.py')
results = json.loads(Path('/scratch/results.json').read_text())
best = results['best_variation']
text = result_path.read_text()
import re
text = re.sub(r'CANONICAL_REPLAY_VARIATION = ".*?"', f'CANONICAL_REPLAY_VARIATION = "{best}"', text, count=1)
result_path.write_text(text)
print(f"Updated canonical replay variation to {best}")
PY
