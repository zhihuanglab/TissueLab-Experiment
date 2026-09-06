set -e
python - <<'PY'
import sys, inspect
sys.path.insert(0, '/shared/lib')
import shared_analysis.stats as stats
import shared_analysis.artifacts as art
print('residualized_loo_predictive_correlation sig:', inspect.signature(stats.residualized_loo_predictive_correlation))
print(inspect.getsource(stats.residualized_loo_predictive_correlation)[:5000])
print('\nartifacts funcs:')
for name in ['write_donor_feature_table','build_results_payload','write_results_payload','validate_results_payload','coerce_results_payload']:
    fn=getattr(art,name)
    print('\n',name, inspect.signature(fn))
    try:
        print(inspect.getsource(fn)[:3500])
    except Exception as e:
        print('src err', e)
PY
