set -e
find /data/autoresearch_runs -path '*round_0001_worker/results.json' | head -n 3 | while read f; do echo 'FILE' $f; cat "$f"; done
