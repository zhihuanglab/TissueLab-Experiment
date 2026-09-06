python - <<'PY'
import os, sys
print('PYTHONPATH env:', os.environ.get('PYTHONPATH'))
print('sys.path head:', sys.path[:8])
import os
print('shared/lib contents:', os.listdir('/shared/lib')[:20])
PY
