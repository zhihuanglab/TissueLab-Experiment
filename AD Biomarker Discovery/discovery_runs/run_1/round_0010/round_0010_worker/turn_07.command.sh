set -euo pipefail
python - <<'PY'
import sys, pandas as pd, numpy as np
sys.path.append('/shared/lib')
from shared_analysis import build_cell_table
slide='/data/H19.33.004-A12-LFB.svs.zarr'
cells=build_cell_table(slide, include_regions=True, include_geometry=True)
print(cells.head().to_string())
print('\ncolumns:', cells.columns.tolist())
print('\nregion counts top:\n', cells['region'].value_counts(dropna=False).head(12).to_string())
print('\ncell types:\n', cells['cell_type'].value_counts().to_string())
ec=cells[cells['region']=='EC']
print('\nEC types:\n', ec['cell_type'].value_counts().to_string())
print('\nEC pyramidal count', ((ec['cell_type']=='Pyramidal Neuron')).sum())
print('EC reactive astro count', ((ec['cell_type']=='Reactive Astrocyte')).sum())
print('area summary EC pyramidal', ec.loc[ec['cell_type']=='Pyramidal Neuron','area'].describe().to_string())
PY
