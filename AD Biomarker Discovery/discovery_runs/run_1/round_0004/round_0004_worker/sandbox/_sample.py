from pathlib import Path
from shared_analysis import load_training_cohort, build_cell_table
cohort = load_training_cohort('/data')
slide = cohort.iloc[0]['slide_name']
cells = build_cell_table(Path('/data')/slide, include_regions=True, include_geometry=False)
print('slide', slide, 'shape', cells.shape)
print('columns', cells.columns.tolist())
print(cells[['cell_type','region','x_centroid','y_centroid']].head().to_string())
print('region counts top:')
print(cells['region'].value_counts(dropna=False).head(12).to_string())
print('cell_type counts top:')
print(cells['cell_type'].value_counts().head(12).to_string())
print('CA1 types top:')
ca1 = cells[cells['region']=='CA1']
print(ca1['cell_type'].value_counts().head(12).to_string())
