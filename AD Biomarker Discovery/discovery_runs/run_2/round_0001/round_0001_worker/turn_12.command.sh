python - <<'PY'
import sys, inspect, textwrap
sys.path.insert(0, '/shared/lib')
import shared_analysis.stats as stats
for name in ['_design_matrix','residualize','_clean_frame']:
    fn=getattr(stats,name)
    print(f'\n--- {name} ---')
    print(textwrap.dedent(inspect.getsource(fn))[:2400])
PY
