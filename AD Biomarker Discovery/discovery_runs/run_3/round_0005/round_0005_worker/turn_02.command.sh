set -e
python - <<'PY'
import json
from pathlib import Path
p=Path('/scratch/context_bundle.json')
obj=json.loads(p.read_text())
print("TOP-LEVEL KEYS:", list(obj.keys()))
for k in ['worker_brief','round','candidate','hypothesis','accepted_panel','current_panel','panel','brief']:
    if k in obj:
        print(f"\n## {k}")
        print(json.dumps(obj[k], indent=2)[:12000])
print("\nSearch keys containing brief/hypothesis/candidate/panel:")
def walk(prefix, x):
    if isinstance(x, dict):
        for k,v in x.items():
            path=f"{prefix}.{k}" if prefix else k
            if any(t in k.lower() for t in ['brief','hypothesis','candidate','panel','round']):
                print(path)
            walk(path,v)
    elif isinstance(x, list):
        for i,v in enumerate(x[:50]):
            walk(f"{prefix}[{i}]",v)
walk("",obj)
PY
if [ -f /scratch/worker_brief.json ]; then
  echo '--- worker_brief.json ---'
  python - <<'PY'
import json
print(json.dumps(json.load(open('/scratch/worker_brief.json')), indent=2)[:20000])
PY
fi
