set -euo pipefail
echo '--- remainder stats helper ---'
sed -n '260,420p' /shared/lib/shared_analysis/stats.py
echo '--- artifacts helper ---'
sed -n '1,260p' /shared/lib/shared_analysis/artifacts.py
echo '--- openslide availability ---'
python - <<'PY'
import importlib.util
print('openslide', importlib.util.find_spec('openslide'))
print('tifffile', importlib.util.find_spec('tifffile'))
PY
