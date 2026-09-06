set -euo pipefail
grep -R "ca1_pyramidal_immune_admixed_reactive_cuff_fraction" -n /data/autoresearch_runs/*/round_0008/round_0008_worker/results.json | head -10
