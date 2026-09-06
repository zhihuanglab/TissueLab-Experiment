set -euo pipefail
python - <<'PY'
from pathlib import Path
path = Path('/scratch/result.py')
text = path.read_text()
text = text.replace(
"        conf_pred_outcome = float(x_test_conf @ beta_outcome)\n        conf_pred_feature = float(x_test_conf @ beta_feature)\n",
"        conf_pred_outcome = float((x_test_conf @ beta_outcome)[0])\n        conf_pred_feature = float((x_test_conf @ beta_feature)[0])\n"
)
text = text.replace(
"        random_state=variation_name.__hash__() & 0xFFFFFFFF,\n",
"        random_state={\"candidate_variant_a\": 11, \"candidate_variant_b\": 17}[variation_name],\n"
)
text = text.replace(
"        random_state=(variation_name.__hash__() + 7) & 0xFFFFFFFF,\n",
"        random_state={\"candidate_variant_a\": 19, \"candidate_variant_b\": 23}[variation_name],\n"
)
text = text.replace(
"""    ranked = sorted(
        variation_metrics,
        key=lambda item: (
            -np.inf if not math.isfinite(float(item["selection_score"])) else -float(item["selection_score"]),
            -np.inf if not math.isfinite(float(item["adjusted_score"])) else -float(item["adjusted_score"]),
        ),
    )
""",
"""    ranked = sorted(
        variation_metrics,
        key=lambda item: (
            not math.isfinite(float(item["selection_score"])),
            -float(item["selection_score"]) if math.isfinite(float(item["selection_score"])) else math.inf,
            -float(item["adjusted_score"]) if math.isfinite(float(item["adjusted_score"])) else math.inf,
        ),
    )
"""
)
path.write_text(text)
PY
python /scratch/result.py
