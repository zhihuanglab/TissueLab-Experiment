#!/usr/bin/env python3
"""Generate source_data/extfig3c_xray_auc_before_after.csv (Extended Data Fig. 3c, NIH ChestX-ray14, 14 conditions).

Run from this folder: python3 make_extfig3c_source_data.py  ->  ../source_data/extfig3c_xray_auc_before_after.csv

Inputs (this folder): before_co-evolving.npy (25,596 x 14 probabilities, NIH condition order, rows in the order of
Q2.5-<Disease>-round1.csv), 'after co-evolving.csv' (per-image probabilities of several classifier heads concatenated, with
image_name), Q2.5-<Disease>-round1.csv (ground_truth per image; slide_name = image id, patient = first 8 characters).

--block first  (default) the FIRST column carrying each condition name in 'after co-evolving.csv' (pandas' plain
               df['Atelectasis']); this is what the published panel and main text use (After Mean AUC 0.8284). For
               Atelectasis, Cardiomegaly, Pneumonia, Pneumothorax, Consolidation and Edema this is the CheXpert-head column.
--block nih    the NIH-head column (third occurrence) for those six conditions; reproduces the 2026-08-29 table exactly.

CI convention (reproduces the 2026-08-29 table to 4 dp): one numpy.random.default_rng(0) stream for the whole table,
rows in plot order (condition in NIH order, before then after); for each row 1,000 image-level resamples
(ci95_case) followed by 1,000 patient-level cluster resamples (ci95_pat, patients sorted by id), percentile 2.5/97.5.
"""
import argparse, os
import numpy as np, pandas as pd
from scipy.stats import rankdata
from sklearn.metrics import roc_auc_score

HERE = os.path.dirname(os.path.abspath(__file__))
NIH = ['Atelectasis', 'Cardiomegaly', 'Effusion', 'Infiltration', 'Mass', 'Nodule', 'Pneumonia', 'Pneumothorax',
       'Consolidation', 'Edema', 'Emphysema', 'Fibrosis', 'Pleural_Thickening', 'Hernia']
B = 1000

def auc(y, s):   # rank AUC with tie averaging = sklearn roc_auc_score for binary labels
    r = rankdata(s); n1 = y.sum(); n0 = len(y) - n1
    return (r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)

def case_ci(rng, y, s):
    n = len(y); out = []
    for _ in range(B):
        i = rng.integers(0, n, n); out.append(auc(y[i], s[i]))
    return np.percentile(out, [2.5, 97.5])

def patient_ci(rng, y, s, pat):
    g, gi = np.unique(pat, return_inverse=True); mem = [np.flatnonzero(gi == k) for k in range(len(g))]; G = len(mem); out = []
    for _ in range(B):
        d = rng.integers(0, G, G); i = np.concatenate([mem[k] for k in d]); out.append(auc(y[i], s[i]))
    return np.percentile(out, [2.5, 97.5])

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--xray-dir', default=HERE)
    ap.add_argument('--block', default='first', choices=('first', 'nih'))
    ap.add_argument('--out', default=os.path.join(HERE, '..', 'source_data', 'extfig3c_xray_auc_before_after.csv'))
    a = ap.parse_args()
    after = pd.read_csv(os.path.join(a.xray_dir, 'after co-evolving.csv'))
    after['slide_name'] = after['image_name'].str.split('/').str[-1].str.replace('.png', '', regex=False)
    before = np.load(os.path.join(a.xray_dir, 'before_co-evolving.npy')).astype(float)
    cols = list(after.columns); rng = np.random.default_rng(0); rows = []
    for i, d in enumerate(NIH):
        gt = pd.read_csv(os.path.join(a.xray_dir, f'Q2.5-{d}-round1.csv'))
        y = gt['ground_truth'].astype(str).str.lower().map({'true': 1, 'false': 0}).values.astype(int)
        pat = gt['slide_name'].astype(str).str[:8].values
        m = gt[['slide_name']].merge(after, on='slide_name', how='left')          # after table in ground-truth row order
        occ = [c for c in cols if c == d or c.startswith(d + '.')]
        col = occ[0] if a.block == 'first' else (occ[2] if len(occ) > 2 else occ[0])
        for system, s in (('Before co-evolution', before[:, i]), ('After co-evolution', m[col].values.astype(float))):
            v = roc_auc_score(y, s); cl, ch = case_ci(rng, y, s); pl, ph = patient_ci(rng, y, s, pat)
            rows.append(dict(plot_order=len(rows) + 1, figure_panel='Extended Data 3c', experiment='Thoracic disease detection',
                             cohort=f'NIH ChestX-ray14, {d}', system=system, metric=f'AUC ({d})', denominator='full_cohort',
                             n_cases=len(y), n_patients=len(np.unique(pat)), n_answered=len(y), n_runs=1,
                             value=round(v, 4), sd_runs='n.a. one run', ci95_case_lo=round(cl, 4), ci95_case_hi=round(ch, 4),
                             ci95_pat_lo=round(pl, 4), ci95_pat_hi=round(ph, 4)))
            print(f'{d:18s} {system:20s} {v:.4f}  [{cl:.4f}, {ch:.4f}]  patients [{pl:.4f}, {ph:.4f}]  column={col}')
    df = pd.DataFrame(rows); os.makedirs(os.path.dirname(a.out), exist_ok=True); df.to_csv(a.out, index=False)
    print(f"mean before {df[df.system.str.startswith('Before')].value.mean():.4f}   mean after {df[df.system.str.startswith('After')].value.mean():.4f}")
    print('wrote', a.out)

if __name__ == '__main__':
    main()
