python - <<'PY'
from pathlib import Path
slide=Path('/data/H19.33.004-A12-LFB.svs')
try:
    import openslide
    s=openslide.OpenSlide(str(slide))
    print('openslide ok')
    for k in ['aperio.MPP','openslide.mpp-x','openslide.mpp-y','aperio.AppMag']:
        print(k, s.properties.get(k))
except Exception as e:
    print('openslide failed', e)
    try:
        import tifffile
        with tifffile.TiffFile(str(slide)) as tf:
            print('tifffile ok pages', len(tf.pages))
            print(tf.pages[0].description[:1000])
    except Exception as e2:
        print('tifffile failed', e2)
PY
