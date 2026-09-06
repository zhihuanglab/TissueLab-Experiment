set -euo pipefail
grep -n "def _adjusted_metrics\|def _selection_score" -n /shared/templates/worker_embedding_mechanistic_template.py
sed -n '80,150p' /shared/templates/worker_embedding_mechanistic_template.py
