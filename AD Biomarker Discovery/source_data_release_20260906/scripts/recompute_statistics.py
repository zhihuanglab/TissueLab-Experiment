"""Run the declared SEA-AD inferential analyses without private data paths."""
import argparse
from pathlib import Path
import subprocess
import sys

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,required=True)
    parser.add_argument('--quick',action='store_true')
    args=parser.parse_args()
    here=Path(__file__).resolve().parent
    inputs=here.parent/'inputs'
    out=args.output_dir.resolve();out.mkdir(parents=True,exist_ok=True)
    perms='200' if args.quick else '100000'
    boots='100' if args.quick else '10000'
    commands=[
        ['run_panel_trajectory_confirmatory_permutation.py','--input-dir',str(inputs),
         '--output',str(out/'panel_permutation.json'),'--permutations',perms],
        ['bootstrap_panel_trajectory_association.py','--input-dir',str(inputs),
         '--output',str(out/'panel_bootstrap.json'),'--bootstrap-samples',boots],
        ['analyze_ad_april_original_statistics.py','--input-dir',str(inputs),
         '--helper-script',str(here/'analyze_ad_reactive_niche_statistics.py'),
         '--output-dir',str(out/'family'),'--n-permutations',perms,'--n-bootstrap',boots,
         '--measurement-stack','NuClass April-compatible reconstruction',
         '--provenance-note','Reconstructed donor tables; see package README for historical/plotting version differences.'],
        ['analyze_ad_method_comparisons.py','--baseline-dir',str(inputs/'baselines'),
         '--feature-matrix',str(inputs/'run_20260420_162449_64fb7a_test.csv'),
         '--output-dir',str(out/'paired_baselines'),'--bootstrap-repetitions',boots,
         '--permutation-repetitions',perms],
    ]
    for command in commands:
        subprocess.run([sys.executable,str(here/command[0]),*command[1:]],check=True)

if __name__=='__main__':main()
