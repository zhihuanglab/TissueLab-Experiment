python - <<'PY'
from pathlib import Path
p = Path('/scratch/result.py')
txt = p.read_text()
txt = txt.replace('import json\nimport math\nimport re\nimport sys\n', 'import json\nimport math\nimport re\nimport sys\nimport warnings\n')
txt = txt.replace('HYPOTHESIS_FAMILY = "CA1 peripyramidal reactive astrocyte lymphocyte-contact area burden"\nCANONICAL_DEFAULT_VARIATION = "candidate_variant_a"\n', 'HYPOTHESIS_FAMILY = "CA1 peripyramidal reactive astrocyte lymphocyte-contact area burden"\nCANONICAL_DEFAULT_VARIATION = "candidate_variant_a"\nFEATURE_NAME = "ca1_peripyramidal_ra_lymphocyte_area_burden_20um"\nFEATURE_COLUMN = "ca1_peripyramidal_ra_lymphocyte_area_burden_20um"\n')
txt = txt.replace('# Make shared helpers importable in this runtime.\nsys.path.insert(0, "/shared/lib")\n', '# Make shared helpers importable in this runtime.\nsys.path.insert(0, "/shared/lib")\nwarnings.filterwarnings("ignore", message="Object at .* is not recognized as a component of a Zarr hierarchy.")\nwarnings.filterwarnings("ignore", message="invalid value encountered in divide")\n')
p.write_text(txt)
PY
python /scratch/result.py >/scratch/last_stdout.txt
tail -n 20 /scratch/last_stdout.txt
