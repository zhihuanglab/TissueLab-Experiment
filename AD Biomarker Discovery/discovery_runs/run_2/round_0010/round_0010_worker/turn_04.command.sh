set -e
echo '--- /data top level ---'
ls -1 /data | sed -n '1,80p'
echo
echo '--- find prior worker result.py under /data/autoresearch_runs ---'
find /data/autoresearch_runs -name result.py | sed -n '1,40p'
