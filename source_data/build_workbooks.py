"""Build one Source Data workbook per figure (one sheet per panel) from the per-panel CSVs in this folder.
Run from source_data/:  python3 build_workbooks.py   ->  workbooks/Source_Data_<figure>.xlsx + workbooks/README.md"""
import glob, os, re
import pandas as pd
here = os.path.dirname(os.path.abspath(__file__)); os.chdir(here)
os.makedirs('workbooks', exist_ok=True)
KIND = {'Figure': ('Fig', 'Fig. '), 'Extended Data': ('ED_Fig', 'ED Fig. '), 'Supplementary Fig.': ('Supp_Fig', 'Supp. Fig. ')}
def fig_key(panel):
    m = re.match(r'(Figure|Extended Data|Supplementary Fig\.) (\d+)([a-z,]*)', panel)
    key, pre = KIND[m.group(1)]
    return f"{key}{m.group(2)}", f"{pre}{m.group(2)}{m.group(3)}"
books, index = {}, []
for f in sorted(glob.glob('*.csv')):
    df = pd.read_csv(f)
    if 'figure_panel' not in df.columns:
        print('skip (no figure_panel column):', f); continue
    for panel, sub in df.groupby('figure_panel', sort=False):
        key, sheet = fig_key(panel)
        books.setdefault(key, []).append((sheet, sub, f))
for key, sheets in books.items():
    out = f'workbooks/Source_Data_{key}.xlsx'
    with pd.ExcelWriter(out, engine='xlsxwriter') as xw:
        for sheet, sub, src in sheets:
            name = sheet[:31]
            sub.to_excel(xw, sheet_name=name, index=False)
            ws = xw.sheets[name]
            for i, c in enumerate(sub.columns):
                w = max([len(str(c))] + [len(str(v)) for v in sub[c].dropna()]) if len(sub) else len(str(c))
                ws.set_column(i, i, max(10, min(60, w + 2)))
            index.append((os.path.basename(out), sheet, src, len(sub)))
    print(f'{out}: {len(sheets)} sheets')
with open('workbooks/README.md', 'w') as fh:
    fh.write('# Source Data workbooks\n\nOne workbook per figure, one sheet per figure panel, as the Nature manuscript checklist requires. '
             'Each sheet is the corresponding per-panel CSV in the parent folder, unchanged. Regenerate with `python3 build_workbooks.py` from `source_data/`.\n\n'
             '| Workbook | Sheet | Source CSV | Rows |\n|---|---|---|---|\n')
    for w, s, src, n in index: fh.write(f'| {w} | {s} | {src} | {n} |\n')
print('total sheets:', len(index))
