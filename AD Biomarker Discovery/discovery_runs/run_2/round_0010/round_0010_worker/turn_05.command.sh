set -e
find /data/autoresearch_runs -path '*run_20260420_215745_90b105*' -name result.py | sed -n '1,50p'
