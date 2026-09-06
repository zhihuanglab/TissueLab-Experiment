python - <<'PY'
from pathlib import Path
path = Path('/scratch/result.py')
text = path.read_text()
text = text.replace(
"    boot_metrics = bootstrap_partial_correlation(\n        feature_table,\n        feature_col=feature_col,\n        outcome_col=OUTCOME_COLUMN,\n        confounds=CONFOUND_COLUMNS,\n        n_boot=2000,\n        random_state=variation_name.__hash__() & 0xFFFFFFFF,\n    )\n    stability_metrics = bootstrap_partial_correlation_stability(\n        feature_table,\n        feature_col=feature_col,\n        outcome_col=OUTCOME_COLUMN,\n        confounds=CONFOUND_COLUMNS,\n        n_boot=400,\n        sample_frac=0.8,\n        random_state=(variation_name.__hash__() + 17) & 0xFFFFFFFF,\n    )\n",
"    seed_base = {\"candidate_variant_a\": 101, \"candidate_variant_b\": 202}[variation_name]\n    boot_metrics = bootstrap_partial_correlation(\n        feature_table,\n        feature_col=feature_col,\n        outcome_col=OUTCOME_COLUMN,\n        confounds=CONFOUND_COLUMNS,\n        n_boot=2000,\n        random_state=seed_base,\n    )\n    stability_metrics = bootstrap_partial_correlation_stability(\n        feature_table,\n        feature_col=feature_col,\n        outcome_col=OUTCOME_COLUMN,\n        confounds=CONFOUND_COLUMNS,\n        n_boot=400,\n        sample_frac=0.8,\n        random_state=seed_base + 17,\n    )\n")
path.write_text(text)
PY
python /scratch/result.py > /scratch/run_stdout.txt
cat /scratch/run_stdout.txt
