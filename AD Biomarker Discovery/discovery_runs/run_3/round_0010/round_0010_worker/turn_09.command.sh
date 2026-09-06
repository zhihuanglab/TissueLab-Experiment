set -euo pipefail
python - <<'PY'
from pathlib import Path
targets = ['ca1_conditional_immune_cuff_fraction_r50um_ra2_ly1',
           'ca1_pyramidal_immune_admixed_reactive_cuff_fraction_r50um_ra2_ly1']
for t in targets:
    print('\nSEARCH', t)
    hits=0
    for p in Path('/data/autoresearch_runs').rglob('results.json'):
        txt = p.read_text(errors='ignore')
        if t in txt:
            print('HIT', p)
            print(txt[:2000])
            hits+=1
            if hits>=3: break
    if not hits:
        print('no hits')
PY
