set -e
find /shared -maxdepth 3 -type d | sed -n '1,40p'
echo
find /shared -maxdepth 4 -type f | sed -n '1,60p'
