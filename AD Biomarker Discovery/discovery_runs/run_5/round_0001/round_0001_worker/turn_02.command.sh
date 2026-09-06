set -euo pipefail
python - <<'PY'
import json, pathlib, textwrap
obj=json.loads(pathlib.Path('/scratch/context_bundle.json').read_text())
print("keys:", list(obj))
for k in ['worker_brief','accepted_panel','current_round','round','candidate','hypothesis','plan']:
    if k in obj:
        print(f"\n--- {k} ---")
        print(json.dumps(obj[k], indent=2)[:12000])
PY
