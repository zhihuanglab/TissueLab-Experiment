set -e
export PYTHONPATH=/shared/lib:${PYTHONPATH}
python - <<'PY'
from pathlib import Path
import inspect
from shared_analysis import build_cell_table, load_training_cohort
import shared_analysis.stats as st
cohort = load_training_cohort('/data')
print("cohort cols:", list(cohort.columns))
slide = Path('/data') / cohort.iloc[0]['slide_name']
print("using slide:", slide.name)
cells = build_cell_table(slide, include_regions=True, include_geometry=False)
print("cell columns:", list(cells.columns))
print("n cells:", len(cells))
for col in ['cell_type','region_label','centroid_x','centroid_y','class_label','label_name','region']:
    if col in cells.columns:
        print(f"\n{col} unique sample:", cells[col].astype(str).value_counts().head(10).to_dict())
print("\npartial_correlation signature:", inspect.signature(st.partial_correlation))
print("leave_one_out_summary signature:", inspect.signature(st.leave_one_out_summary))
print("residualized_loo_predictive_correlation signature:", inspect.signature(st.residualized_loo_predictive_correlation))
PY
