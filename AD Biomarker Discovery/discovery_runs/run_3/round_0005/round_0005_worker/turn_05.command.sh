set -e
echo 'Shared lib contents:'
find /shared -maxdepth 3 -type f | sed -n '1,80p'
echo 'Python path check:'
python - <<'PY'
import sys
print(sys.path)
PY
