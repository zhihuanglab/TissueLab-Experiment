set -euo pipefail
python - <<'PY'
import json, textwrap, os, pathlib
for path in ["/scratch/context_bundle.json", "/scratch/program.md", "/scratch/result.py"]:
    print("\n" + "="*20, path, "="*20)
    if path.endswith(".json"):
        with open(path) as f:
            obj = json.load(f)
        print(json.dumps(obj, indent=2)[:30000])
    else:
        with open(path) as f:
            txt = f.read()
        print(txt[:30000])
PY
ls -lah /data | sed -n '1,120p'
