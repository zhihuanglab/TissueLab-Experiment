set -euo pipefail
echo 'shared_analysis package files:'
find /shared/lib/shared_analysis -maxdepth 2 -type f | sed -n '1,80p'
echo '--- __init__.py ---'
sed -n '1,200p' /shared/lib/shared_analysis/__init__.py
echo '--- maybe cell table implementation file names ---'
grep -R "def build_cell_table" -n /shared/lib/shared_analysis
