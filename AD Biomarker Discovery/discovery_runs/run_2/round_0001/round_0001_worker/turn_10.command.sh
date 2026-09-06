python - <<'PY'
import sys, inspect, textwrap, re
sys.path.insert(0, '/shared/lib')
import shared_analysis.sea_ad_lfb as mod
src = inspect.getsource(mod)
for pat in ['mpp','openslide','tiff','slide']:
    if pat in src.lower():
        pass
# Print relevant snippets around mpp
for m in re.finditer('mpp', src, flags=re.I):
    start=max(0,m.start()-250); end=min(len(src),m.end()+400)
    print('--- snippet ---')
    print(src[start:end])
    print()
PY
