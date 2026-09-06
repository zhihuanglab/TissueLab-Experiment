set -euo pipefail
python - <<'PY'
try:
    import openslide
    print('openslide version ok')
    s=openslide.OpenSlide('/data/H19.33.004-A12-LFB.svs')
    print(s.properties.get('openslide.mpp-x'), s.properties.get('aperio.MPP'))
except Exception as e:
    print('openslide failed', repr(e))
PY
