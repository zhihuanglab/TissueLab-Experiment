set -e
echo 'shared tree'
find /shared -maxdepth 3 -type d | sed -n '1,60p'
echo 'scratch files'
ls -la /scratch | sed -n '1,80p'
