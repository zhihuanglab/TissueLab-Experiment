set -e
python - <<'PY'
import sys, inspect
sys.path.insert(0, '/shared/lib')
import shared_analysis.stats as stats
for name in ['leave_one_out_partial_correlation','residualize','_design_matrix','bootstrap_partial_correlation']:
    fn=getattr(stats,name)
    print('\n###',name, inspect.signature(fn))
    print(inspect.getsource(fn)[:4000])
PY
