set -euo pipefail
sed -n '260,380p' /shared/lib/shared_analysis/stats.py
echo '--- artifacts.py ---'
sed -n '1,240p' /shared/lib/shared_analysis/artifacts.py
