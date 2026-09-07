"""Reconstruct the seven numerical panels from the public CSVs."""
import argparse
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[3]
DATA=ROOT/'source_data'
COLORS=['#b1c0e1','#7fc6ab']

def finish(fig,ax,title,out,name):
    ax.set_title(title)
    ax.spines[['top','right']].set_visible(False)
    ax.grid(axis='y',alpha=.15);ax.set_axisbelow(True)
    fig.tight_layout();fig.savefig(out/(name+'.png'),dpi=180);plt.close(fig)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output-dir',type=Path,required=True);args=p.parse_args()
    out=args.output_dir;out.mkdir(parents=True,exist_ok=True)
    d=pd.read_csv(DATA/'extfig6b_seaad_panel_trajectory.csv')
    fig,ax=plt.subplots(figsize=(7,4.6))
    test=d[d.cohort.eq('SEA-AD test')]
    ax.fill_between(test['round'],0,test.null_p95,color='#e4e6ee',label='Held-out null 95th percentile')
    ax.plot(test['round'],test.null_median,':',color='grey',label='Held-out null median')
    for color,(cohort,g) in zip(COLORS,d.groupby('cohort',sort=False)):
        ax.plot(g['round'],g.value,'o-',color=color,label=cohort)
        ax.fill_between(g['round'],g.value-g.sd_runs,g.value+g.sd_runs,color=color,alpha=.25)
    ax.set(xlabel='Discovery round',ylabel='Multiple partial R',ylim=(0,.85));ax.legend(fontsize=8)
    finish(fig,ax,'Accumulated panel association',out,'extfig6b')

    d=pd.read_csv(DATA/'extfig6c_seaad_biomarker_comparison.csv')
    fig,ax=plt.subplots(figsize=(10,4.8));x=np.arange(len(d))
    ax.bar(x,d.value,yerr=d.sd_subsample,color=[COLORS[0]]*4+['#e8d8a8']+[COLORS[1]]*3,capsize=3)
    names=['UNI\nABMIL','TITAN\nlinear probe','MUSK\nABMIL','TITAN\nfine-tuning','CA1 GFAP\narea','Reactive\nfraction','Reactive/astroglial\nratio, 35 um','Peripyramidal\nreactivity, 35 um']
    ax.set_xticks(x,names,fontsize=8);ax.axhline(0,color='black',lw=.6);ax.axvline(3.5,ls=':',color='grey')
    ax.set_ylabel('Direction-aligned partial r')
    finish(fig,ax,'Held-out memory-decline association',out,'extfig6c')

    for filename,group,title,name in [
        ('extfig7k_seaad_adni_dementia.csv','cohort','Clinical dementia status','extfig7k'),
        ('suppfig4b_seaad_tool_transfer.csv','measurement_stack','Measurement-stack comparison','suppfig4b')]:
        d=pd.read_csv(DATA/filename);fig,ax=plt.subplots(figsize=(7,4.4));measurements=list(d.system.unique());x=np.arange(3)
        for j,(arm,g) in enumerate(d.groupby(group,sort=False)):
            g=g.set_index('system').loc[measurements]
            ax.bar(x+(j-.5)*.35,g.value,.35,yerr=g.sd_subsample,color=COLORS[j],capsize=3,label=arm)
        ax.set_xticks(x,['Reactive\nfraction','Reactive/astroglial\nratio, 35 um','Peripyramidal\nreactivity, 35 um'])
        ax.set_ylabel('Partial r' if name=='extfig7k' else 'Absolute partial r');ax.legend(fontsize=8)
        finish(fig,ax,title,out,name)

    d=pd.read_csv(DATA/'suppfig4a_classifier_training_labels.csv');fig,ax=plt.subplots(figsize=(6.5,5));bottom=np.zeros(2)
    for cls,g in d.groupby('cell_class',sort=False):
        g=g.set_index('measurement_stack').loc[['NuClass','NeuroFM']]
        ax.bar([0,1],g.label_count,bottom=bottom,color=g.color_hex.iloc[0],label=cls)
        bottom+=g.label_count.to_numpy()
    ax.set_xticks([0,1],['NuClass','NeuroFM']);ax.set_ylabel('Pathologist labels');ax.legend(fontsize=7,bbox_to_anchor=(1.01,1))
    finish(fig,ax,'Classifier training labels',out,'suppfig4a')

    d=pd.read_csv(DATA/'suppfig5a_cognitive_followup.csv');fig,axes=plt.subplots(1,2,figsize=(8,4))
    for ax,(variable,g) in zip(axes,d.groupby('variable',sort=False)):
        stats=[dict(label=r.cohort,q1=r.q1,med=r.median,q3=r.q3,whislo=r.whisker_low,whishi=r.whisker_high,fliers=[]) for r in g.itertuples()]
        ax.bxp(stats,showfliers=False);ax.set_title(variable);ax.set_ylabel(g.unit.iloc[0])
    fig.suptitle('Follow-up summaries (ADNI individual points not redistributed)');fig.tight_layout();fig.savefig(out/'suppfig5a.png',dpi=180);plt.close(fig)

    d=pd.read_csv(DATA/'suppfig5b_adni_followup_associations.csv');fig,ax=plt.subplots(figsize=(8,4.5));groups=['all','>=5y','>=7y','held-out reference'];x=np.arange(4)
    for j,(measurement,g) in enumerate(d.groupby('system',sort=False)):
        g=g.set_index('group').loc[groups];ax.bar(x+(j-1)*.24,g.value,.24,color=['#d9eee8','#a8dcc6','#6cb494'][j],label=measurement)
    ax.set_xticks(x,['ADNI all\nn=56','ADNI >=5 y\nn=25','ADNI >=7 y\nn=17','SEA-AD held-out\nn=34'])
    ax.set_ylabel('Direction-aligned partial r');ax.legend(fontsize=7)
    finish(fig,ax,'Memory slope and follow-up duration',out,'suppfig5b')

if __name__=='__main__':main()
