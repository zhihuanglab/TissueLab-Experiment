set -euo pipefail
export PYTHONPATH=/shared/lib:${PYTHONPATH:-}
python - <<'PY'
import inspect
from shared_analysis import stats, artifacts, sea_ad_lfb
for mod, names in [
    (stats, ['partial_correlation','leave_one_out_summary','residualized_loo_predictive_correlation']),
    (artifacts, ['write_donor_feature_table','build_results_payload','write_results_payload','validate_results_payload']),
    (sea_ad_lfb, ['build_cell_table','load_region_polygons']),
]:
    print(f'## MODULE {mod.__name__}')
    for name in names:
        fn=getattr(mod,name)
        print(name, inspect.signature(fn))
        doc=inspect.getdoc(fn)
        if doc:
            print(' doc:', doc.split('\\n')[0])
    print()
print('--- build_cell_table one slide ---')
from pathlib import Path
from shared_analysis import load_training_cohort, build_cell_table
cohort = load_training_cohort('/data')
slide = cohort.iloc[0]['slide_name']
cells = build_cell_table(Path('/data')/slide, include_regions=True, include_geometry=False)
print('slide', slide, 'shape', cells.shape)
print('columns', cells.columns.tolist())
print(cells[['cell_type','region','x_centroid','y_centroid']].head().to_string())
print('region counts top:')
print(cells['region'].value_counts(dropna=False).head(10).to_string())
print('cell_type counts top:')
print(cells['cell_type'].value_counts().head(10).to_string())
PY
