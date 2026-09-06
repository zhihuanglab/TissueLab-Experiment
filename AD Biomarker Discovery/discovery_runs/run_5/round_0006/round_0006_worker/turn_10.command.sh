python - <<'PY'
from pathlib import Path
path = Path('/scratch/result.py')
text = path.read_text()
text = text.replace(
    '        feat_resid_test = feat_test - float(x_test @ beta_feat)\n'
    '        out_resid_test = out_test - float(x_test @ beta_out)\n\n'
    '        pred_resid = slope * feat_resid_test if np.isfinite(slope) else float("nan")\n'
    '        pred_raw = float(x_test @ beta_out) + pred_resid if np.isfinite(pred_resid) else float("nan")\n',
    '        feat_resid_test = feat_test - float((x_test @ beta_feat)[0])\n'
    '        out_resid_test = out_test - float((x_test @ beta_out)[0])\n\n'
    '        pred_resid = slope * feat_resid_test if np.isfinite(slope) else float("nan")\n'
    '        pred_raw = float((x_test @ beta_out)[0]) + pred_resid if np.isfinite(pred_resid) else float("nan")\n'
)
text = text.replace(
    '    ranked_variations.sort(\n'
    '        key=lambda d: (\n'
    '            -np.inf if not np.isfinite(float(d["selection_score"])) else -float(d["selection_score"]),\n'
    '            -np.inf if not np.isfinite(float(d["adjusted_score"])) else -float(d["adjusted_score"]),\n'
    '            -np.inf if not np.isfinite(float(d["loo_predictive_r"])) else -abs(float(d["loo_predictive_r"])),\n'
    '        )\n'
    '    )\n',
    '    ranked_variations.sort(\n'
    '        key=lambda d: (\n'
    '            float("-inf") if not np.isfinite(float(d["selection_score"])) else float(d["selection_score"]),\n'
    '            float("-inf") if not np.isfinite(float(d["adjusted_score"])) else float(d["adjusted_score"]),\n'
    '            float("-inf") if not np.isfinite(float(d["loo_predictive_r"])) else abs(float(d["loo_predictive_r"])),\n'
    '        ),\n'
    '        reverse=True,\n'
    '    )\n'
)
path.write_text(text)
PY
python /scratch/result.py
echo '--- results.json ---'
python - <<'PY'
import json
with open('/scratch/results.json') as f:
    data=json.load(f)
print(json.dumps(data, indent=2)[:2000])
PY
echo '--- report.md ---'
sed -n '1,240p' /scratch/report.md
