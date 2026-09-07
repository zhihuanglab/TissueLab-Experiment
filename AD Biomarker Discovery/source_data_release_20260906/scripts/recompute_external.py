"""Compute aggregates from authorized local ADNI inputs; never export donor rows."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd

def partial(x,y,z):
    rx=x-z@np.linalg.lstsq(z,x,rcond=None)[0]
    ry=y-z@np.linalg.lstsq(z,y,rcond=None)[0]
    return float(np.corrcoef(rx,ry)[0,1])

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input',type=Path,required=True,
        help='Authorized CSV: age_at_death,sex_binary,braak_numeric,cerad_ordinal,span_years,bulk,hit,peri; slope_adni_mem_mixed and/or dementia_status')
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();d=pd.read_csv(args.input)
    covs=['age_at_death','sex_binary','braak_numeric','cerad_ordinal']
    rows=[]
    for outcome in ['slope_adni_mem_mixed','dementia_status']:
        if outcome not in d:continue
        for threshold in ([0,5,7] if outcome=='slope_adni_mem_mixed' else [0]):
            for f in ['bulk','hit','peri']:
                sub=d[d.span_years.ge(threshold)].dropna(subset=covs+[outcome,f])
                if len(sub)<=len(covs)+2:raise ValueError('Insufficient complete cases')
                z=np.column_stack([np.ones(len(sub)),sub[covs].to_numpy(float)])
                r=partial(sub[f].to_numpy(float),sub[outcome].to_numpy(float),z)
                rows.append(dict(outcome=outcome,minimum_followup_years=threshold,
                    measurement=f,n=len(sub),signed_partial_r=r,
                    direction_aligned_partial_r=(-r if outcome=='slope_adni_mem_mixed' else r)))
    if not rows:raise ValueError('No supported outcome column supplied')
    args.output.parent.mkdir(parents=True,exist_ok=True)
    pd.DataFrame(rows).to_csv(args.output,index=False)

if __name__=='__main__':main()
