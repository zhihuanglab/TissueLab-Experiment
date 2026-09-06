set -euo pipefail
echo '--- worker_analysis_template.py ---'
sed -n '1,260p' /shared/templates/worker_analysis_template.py
echo '--- shared_analysis __init__.py ---'
sed -n '1,220p' /shared/lib/shared_analysis/__init__.py
echo '--- shared_analysis stats.py ---'
sed -n '1,260p' /shared/lib/shared_analysis/stats.py
echo '--- round8 result.py ---'
sed -n '1,260p' /data/autoresearch_runs/run_20260421_014751_94b13b/round_0008/round_0008_worker/result.py
echo '--- round8 results.json excerpt ---'
sed -n '1,240p' /data/autoresearch_runs/run_20260421_014751_94b13b/round_0008/round_0008_worker/results.json
