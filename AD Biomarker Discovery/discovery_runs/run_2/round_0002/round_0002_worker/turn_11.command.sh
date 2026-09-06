set -euo pipefail
python - <<'PY'
from pathlib import Path
p = Path('/scratch/result.py')
text = p.read_text()
text = text.replace(
    '        confound_pred_test = float(x_test_conf @ beta_outcome)\n',
    '        confound_pred_test = float((x_test_conf @ beta_outcome)[0])\n',
)
text = text.replace(
    '        resid_feature_test = float(test[feature_col].to_numpy(dtype=float) - x_test_conf @ beta_feature)\n',
    '        resid_feature_test = float((test[feature_col].to_numpy(dtype=float) - x_test_conf @ beta_feature)[0])\n',
)
text = text.replace(
    '    for row in best_loo.itertuples(index=False):\n'
    '        print(\n'
    '            f"  {row.donor_id}  {row.outcome:.4f}  {row.predicted:.4f}  {getattr(row, best_feature_col):.6f}"\n'
    '        )\n',
    '    for _, row in best_loo.iterrows():\n'
    '        print(\n'
    '            f"  {row[\'donor_id\']}  {row[\'outcome\']:.4f}  {row[\'predicted\']:.4f}  {row[best_feature_col]:.6f}"\n'
    '        )\n',
)
p.write_text(text)
PY
python /scratch/result.py
ls -l /scratch/result.py /scratch/results.json /scratch/report.md /scratch/donor_feature_table.csv
