python - <<'PY'
from pathlib import Path
p=Path('/scratch/result.py')
text=p.read_text()
text=text.replace('CANONICAL_VARIATION_NAME = "k2_within_30um"', 'CANONICAL_VARIATION_NAME = "k3_within_30um"')
text=text.replace('FEATURE_NAME = "ca1_pyramidal_reactive_astro_crowding_fraction_k2_30um"', 'FEATURE_NAME = "ca1_pyramidal_reactive_astro_crowding_fraction_k3_30um"')
text=text.replace('FEATURE_COLUMN = "ca1_pyramidal_reactive_astro_crowding_fraction_k2_30um"', 'FEATURE_COLUMN = "ca1_pyramidal_reactive_astro_crowding_fraction_k3_30um"')
# add warning filter if not present
needle='import json\nimport math\nimport sys\nfrom pathlib import Path\n\nsys.path.append("/shared/lib")\n'
repl='import json\nimport math\nimport sys\nimport warnings\nfrom pathlib import Path\n\nwarnings.filterwarnings("ignore", message="Object at .* is not recognized as a component of a Zarr hierarchy.")\nsys.path.append("/shared/lib")\n'
text=text.replace(needle,repl)
p.write_text(text)
PY
python /scratch/result.py > /scratch/run_stdout.txt 2>/scratch/run_stderr.txt
head -n 15 /scratch/run_stdout.txt
echo 'stderr lines:' $(wc -l < /scratch/run_stderr.txt)
