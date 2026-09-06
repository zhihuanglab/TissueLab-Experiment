set -euo pipefail
python - <<'PY'
import sys, inspect, textwrap
sys.path.insert(0, '/shared/lib')
import shared_analysis.stats as st
import shared_analysis.artifacts as art
import shared_analysis.sea_ad_lfb as sea
print('--- stats.leave_one_out_summary ---')
print(textwrap.dedent(inspect.getsource(st.leave_one_out_summary))[:5000])
print('--- stats.partial_correlation ---')
print(textwrap.dedent(inspect.getsource(st.partial_correlation))[:3000])
print('--- artifacts write helpers ---')
for name in ['write_donor_feature_table','build_results_payload','write_results_payload','coerce_results_payload','validate_results_payload']:
    print(f'FUNCTION {name}')
    print(textwrap.dedent(inspect.getsource(getattr(art,name)))[:2500])
    print()
print('--- sea_ad_lfb.build_cell_table ---')
print(textwrap.dedent(inspect.getsource(sea.build_cell_table))[:5000])
PY
