set -euo pipefail
cat > /scratch/report.md <<'MD'
## Summary
Tested CA1 pyramidal immune-admixed severe reactive-cuff fractions at 35 um versus 50 um; the 50 um variation won with selection score 0.4752.

## Metrics
Winning variation: `ca1_pyramidal_immune_admixed_reactive_cuff_fraction_r50um_ra2_ly1`

- IS partial r: -0.4752
- Selection score: 0.4752
- LOO predictive r: 0.2352
- IS-LOO Gap: 0.2400
- Adjusted Score: 0.2352

Other tested variation:

- `ca1_pyramidal_immune_admixed_reactive_cuff_fraction_r35um_ra2_ly1`: partial r 0.1150, selection score 0.1150, LOO predictive r -0.0445

Local ranking therefore favored the broader 50 um immune-admixed cuff strongly over the tighter 35 um version, which was too weak and directionally unstable.

## Findings
1. What worked and why  
   The biologically specific signal came from a CA1 peripyramidal inflammatory niche: the donor-level fraction of CA1 pyramidal neurons with a severe reactive-astrocyte cuff (`>=2` reactive astrocytes) that also included at least one nearby lymphocyte within 50 um. This won because it preserved the previously strong CA1 pyramidal reactive-cuff idea from earlier rounds, but added a distinct inflammatory component that appears to track a subset of donors with worse memory decline.

2. What failed and why  
   The 35 um version failed. Its feature was too sparse: many donors had zero lymphocyte-admixed severe cuffs at that tighter radius, so the signal became near-binary and unstable. That likely explains the weak in-sample association (partial r 0.1150) and negative LOO diagnostic (-0.0445). In short, the immune component exists, but it is not consistently packed tightly enough to be captured at 35 um.

3. Error pattern: which donors are consistently wrong and what they share  
   The largest misses split into two biologically meaningful groups.  
   - Undercalled severe-decline donors with little or no immune-admixed cuff signal: `H20.33.018`, `H20.33.038`, `H20.33.046`, and `H21.33.023` all had more negative memory-decline outcomes than predicted despite very low or zero feature values. They appear to have reactive CA1 pathology without much lymphocyte admixture, so this marker misses non-immune reactive cuffs.  
   - Overcalled inflammatory donors with only moderate decline: `H21.33.043`, `H21.33.044`, and `H20.33.013` had relatively higher immune-admixed cuff fractions and CA1 lymphocyte counts, but less negative outcomes than predicted. They may reflect inflammatory involvement that is real histologically but not tightly coupled to memory decline in the same way across donors.

## Rationale
The best variation is biologically coherent because it keeps the population and niche already implicated by the accepted panel—CA1 pyramidal neurons with reactive astrocytic cuffing—then asks whether the more inflammatory subset of those cuffs carries distinct signal. The 50 um radius beat 35 um because lymphocytes are rarer and probably sit a bit farther from the pyramidal soma than the astrocytes themselves; the broader radius captures that admixed inflammatory microenvironment without changing the underlying niche identity.

This looks potentially additive in interpretation, even if not obviously stronger than the best astrocyte-only cuff features by itself: it isolates an inflammatory subclass of the CA1 reactive cuff rather than re-measuring generic reactive burden.

## Interpretation
The signal seems to mean that memory decline is worse in donors where CA1 pyramidal neurons are embedded in a broader perineuronal inflammatory cuff, not just a purely astrocytic one.

- Population: CA1 Pyramidal Neuron targets, conditioned on nearby Reactive Astrocytes and Lymphocytes
- Niche: CA1 peripyramidal inflammatory cuff within 50 um
- Feature summary: fraction of all CA1 pyramidal neurons with `>=2` nearby Reactive Astrocytes and `>=1` nearby Lymphocyte
- Simplest observable tissue pattern: CA1 neurons sitting inside a dense reactive astrocyte halo that also contains at least one small lymphocyte in the surrounding perisomatic neighborhood

## Next
Stay in the same CA1 pyramidal inflammatory-cuff family, but reduce zero inflation by keeping the 50 um winner and locally sweeping the immune threshold rather than the radius—for example compare `>=1` versus lymphocyte-count-weighted exposure among already severe reactive-cuffed neurons—to test whether the main limitation is sparsity of the hard lymphocyte gate.
MD
ls -l /scratch/result.py /scratch/report.md /scratch/results.json /scratch/donor_feature_table.csv /scratch/best_variation.json
