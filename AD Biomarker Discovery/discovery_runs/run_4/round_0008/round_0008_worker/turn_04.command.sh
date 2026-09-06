set -euo pipefail
export PYTHONPATH=/shared/lib:${PYTHONPATH:-}
python - <<'PY'
from pathlib import Path
import inspect
from shared_analysis import build_cell_table, load_training_cohort
import shared_analysis.sea_ad_lfb as sea
import shared_analysis.stats as st
print('build_cell_table signature:', inspect.signature(build_cell_table))
print('load_region_polygons signature:', inspect.signature(sea.load_region_polygons))
print('partial_correlation signature:', inspect.signature(st.partial_correlation))
print('residualized_loo_predictive_correlation signature:', inspect.signature(st.residualized_loo_predictive_correlation))
cohort=load_training_cohort('/data')
slide_name=cohort.iloc[0]['slide_name']
print('sample slide', slide_name)
cells=build_cell_table(Path('/data')/slide_name, include_regions=True, include_geometry=False)
print('columns:', list(cells.columns))
print('n=', len(cells))
for col in ['label','predicted_label','class_label','region','x_centroid','y_centroid','centroid_x','centroid_y']:
    if col in cells.columns:
        print(col, cells[col].head().tolist()[:5])
label_col = next((c for c in ['label','predicted_label','class_label'] if c in cells.columns), None)
print('label col:', label_col)
print('label value counts head:', cells[label_col].value_counts().head(10).to_dict())
print('region counts:', cells['region'].value_counts().head(10).to_dict())
PY
