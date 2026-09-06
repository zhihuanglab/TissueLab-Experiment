set -euo pipefail
sed -n '1,280p' /shared/lib/shared_analysis/stats.py
echo '---'
sed -n '280,430p' /shared/lib/shared_analysis/artifacts.py
