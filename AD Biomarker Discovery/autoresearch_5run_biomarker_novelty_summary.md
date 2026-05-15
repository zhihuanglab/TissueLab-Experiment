# Autoresearch 5-Run Biomarker Summary

## Scope

This summary covers the **5 10-round autoresearch runs**:

- `run_20260420_162449_64fb7a`
- `run_20260420_215745_90b105`
- `run_20260421_001115_585b70`
- `run_20260421_014751_94b13b`
- `run_20260421_052101_efbee1`

## Aggregate Panel Result

Across the 5 runs, the average panel improved steadily on both train and held-out test:

- round 1 mean train/test panel `r`: `0.373 / 0.422`
- round 5 mean train/test panel `r`: `0.594 / 0.598`
- round 10 mean train/test panel `r`: `0.680 / 0.689`

This is the main take-home result. The system shows a coherent non-collapsing train/test trajectory across independent 10-round searches.

## Consensus Biological Theme

The runs converged on the same broad biological story:

1. **CA1 pyramidal vulnerability**
2. **reactive astrocyte burden / enrichment**
3. **local reactive-astrocyte organization around CA1 pyramidal neurons**
4. **later niche refinements involving immune cells, corpora amylacea, and sometimes oligodendrocyte depletion**

In plain terms, the repeated result was not merely "more pathology in AD tissue." The system kept returning to a more specific hypothesis:

> memory decline is associated with a **CA1-centered reactive astrocyte microenvironment**, especially when reactive astrocytes form local cuffs, crowded niches, or hotspot-like neighborhoods around vulnerable pyramidal neurons.

## What Was Reproducibly Rediscovered

### 1. CA1 pyramidal vulnerability

This is the least novel part, but it was repeatedly recovered:

- `ca1_pyramidal_fraction`
- `ca1_pyramidal_burden__ca1_pyramidal_fraction`
- `ca1_pyramidal_near_reactive_astro_fraction_30um`
- multiple later derivatives based on pyramidal exclusion, depletion, or atrophy near reactive astrocytes

Interpretation:

- This is a real result, but it is fundamentally confirmatory.
- It matches long-established literature showing selective CA1 / subicular pyramidal vulnerability in AD.

Representative literature:

- Davies et al., 1992: [The effect of age and Alzheimer's disease on pyramidal neuron density in the individual fields of the hippocampal formation](https://pubmed.ncbi.nlm.nih.gov/1621507/)

Publication value:

- **Publication-worthy as confirmation and prioritization**
- **Not publication-worthy as a novelty claim by itself**

### 2. Reactive astrocyte burden and enrichment

This was the most stable recurring feature class across runs:

- `ca1_reactive_astrocyte_lineage_fraction`
- `ca1_reactive_astrocyte_enrichment`
- `astrocytic_reactivity_fraction`
- `ca1_reactive_astro_fraction`

These features were among the strongest standalone discoveries:

- repeated `partial_r` around `-0.49`
- repeated nominal `p ~ 0.0027`
- full coverage in all clean runs

Interpretation:

- Reactive astrogliosis is a genuine and replicated signal in this dataset.
- This is also not novel in the disease-mechanism sense, but it is one of the strongest pillars of the discovered panels.

Representative literature:

- Vijayan et al., 1991: [Astrocyte hypertrophy in the Alzheimer's disease hippocampal formation](https://pubmed.ncbi.nlm.nih.gov/2013308/)
- Rodríguez-Arellano et al., 2016: [Astrocytes in physiological aging and Alzheimer's disease](https://pubmed.ncbi.nlm.nih.gov/25595973/)

Publication value:

- **Definitely publication-worthy as a reproducible panel component**
- **Not a novel standalone biological claim**

### 3. Peripyramidal reactive-astrocyte neighborhoods, cuffs, and crowding

This is where the runs become more interesting.

Recurring examples:

- `ca1_pyramidal_reactive_astro_niche_fraction_80px`
- `ca1_pyramidal_reactive_astro_coverage_60px`
- `ca1_pyramidal_reactive_niche__reactive_over_astroglial_neighbored_pyramidal_fraction_r70px`
- `peripyramidal_reactivity_r35um`
- `ca1_pyramidal_cuffed_fraction_r35um_ge1`
- `ca1_pyramidal_cuffed_fraction_r35um_ge2`
- `ca1_pyramidal_supermajority_reactive_cuff_fraction_r35um_ge2`
- `ca1_pyramidal_reactive_astro_crowding_fraction_k3_30um`
- `ca1_pyramidal_reactive_astro_isolated_fraction_k5_30um`

These were not random one-offs. They recurred independently across runs with strong standalone support in several cases:

- niche fraction / coverage terms in the early runs had `p ~ 0.001-0.003`
- cuffing/crowding terms in the overnight runs also had `p ~ 0.001-0.003`

Interpretation:

- This is the strongest scientifically interesting refinement produced by the system.
- The repeated pattern is not just "reactive astrocytes are abundant."
- It is specifically that **reactive astrocytes become organized around CA1 pyramidal neurons in a way that tracks memory decline**.

Novelty judgment:

- The component biology is plausible and connected to known reactive astrocyte pathology.
- But the exact **donor-level biomarker formulation around peripyramidal reactive cuffs / crowding / neighborhood fractions** looks substantially less standard than plain astrocyte burden.

Publication value:

- **Publication-worthy**
- **One of the genuinely novel-ish outputs of these runs**

## More Tentative but Interesting Findings

### 4. Reactive-astro niche morphology and hypertrophy

Recurring examples:

- `ca1_reactive_astrocyte_hypertrophy`
- `ca1_reactive_astrocyte_hypertrophy__reactive_minus_astrocyte_median_log1p_area_ca1`
- `ca1_peripyramidal_reactive_astro_area_median_30um_um2`
- `ca1_peripyramidal_reactive_area_median_r35um`
- proximal small-pyramidal / lower-tail area terms in `run_20260420_162449_64fb7a`

These had mixed standalone strength:

- some hypertrophy/area terms were nominally significant (`p ~ 0.026-0.035`)
- some later morphology refinements were weak alone but still improved the panel

Interpretation:

- Astrocyte morphology and local neuronal size-shift are likely meaningful refinements of the core signal.
- They are interesting, but less replicated cleanly than the burden/enrichment and cuffing/crowding classes.

Publication value:

- **Worth reporting as supporting panel structure**
- **Best framed as refinement, not the headline finding**

### 5. Immune / lymphocyte-associated reactive niches

Recurring examples:

- `ca1_pyramidal_immune_admixed_reactive_cuff_fraction_r50um_ra2_ly1`
- `ca1_conditional_immune_cuff_fraction_r50um_ra2_ly1`
- `ca1_conditional_immune_cuff_fraction_r50um_ra2_ly2`
- `ca1-reactive-lymphocyte-niche`
- `ca1_peripyramidal_reactive_astro_lymphocyte_contact_fraction_20um`
- `ca1_pyramidal_reactive_astro_lymphocyte_triad_fraction_25um`
- `ca1_ra_exposed_pyramidal_lym_ge1_fraction_25um`
- `ca1_peripyramidal_ra_lymphocyte_area_burden_20um`
- the hotspot-local lymphocyte terms in `run_20260420_215745_90b105`

This class is important, but the statistical pattern is different:

- some immune-linked terms were individually significant:
  - `ca1_pyramidal_immune_admixed_reactive_cuff_fraction_r50um_ra2_ly1`: `p = 0.0039`
  - `ca1-reactive-lymphocyte-niche`: `p = 0.0398`
- many later lymphocyte terms were **not** individually strong:
  - several had `p > 0.2`
  - several had modest standalone `partial_r`

However, they still improved panel performance.

Interpretation:

- These behave more like **panel-complement features** than universal standalone biomarkers.
- Scientifically, this is still interesting because the system repeatedly converged on immune context around reactive astrocyte niches.

Representative literature:

- Zeng et al., 2024: [T cell infiltration mediates neurodegeneration and cognitive decline in Alzheimer's disease](https://pubmed.ncbi.nlm.nih.gov/38437992/)

Novelty judgment:

- Adaptive immune involvement in AD is not novel.
- What is more unusual here is the repeated emergence of **lymphocyte-enriched reactive-astrocyte niches specifically in CA1 panel construction**.

Publication value:

- **Publication-worthy as an exploratory panel-level result**
- **Too early to frame as a validated standalone biomarker class**

### 6. Corpora amylacea-associated reactive niches

Recurring examples:

- `ca1_reactive_hotspot_corpora_amylacea_enrichment__reactive_minus_astro_hotspot_local_ca_fraction_r70px_k3`
- `ca1_reactive_corpora_amylacea_neighbor_fraction_r96px`

These were less replicated and individually weaker:

- one corpora-amylacea feature had `p = 0.136`
- another had `p = 0.490`

Interpretation:

- These are interesting because they suggest that the waste-handling / degeneration-associated compartment might modulate the reactive niche signal.
- But they are not yet strong enough to be treated as a core confirmed result.

Representative literature:

- Kang et al., 2022: [Corpora amylacea are associated with tau burden and cognitive status in Alzheimer's disease](https://pubmed.ncbi.nlm.nih.gov/35941704/)
- Dallmeier et al., 2024: [Corpora amylacea negatively correlate with hippocampal tau pathology in Alzheimer's disease](https://pubmed.ncbi.nlm.nih.gov/38486969/)

Publication value:

- **Interesting exploratory subfinding**
- **Not strong enough to headline**

### 7. Oligodendrocyte gap / depletion around reactive niches

Examples:

- `ca1_reactive_hotspot_oligodendrocyte_depletion`
- `ca1_reactive_oligodendrocyte_gap_r80px`

This class is mixed:

- one oligodendrocyte feature was individually significant (`p = 0.018`)
- another was weak (`p = 0.439`)

Interpretation:

- Oligodendroglial involvement in AD is plausible and increasingly supported in the literature.
- In these runs, oligodendrocyte-linked features are interesting but do not recur as strongly as the astrocyte / pyramidal / immune terms.

Representative literature:

- Butt et al., 2019: [Oligodendroglial Cells in Alzheimer's Disease](https://pubmed.ncbi.nlm.nih.gov/31583593/)
- Clayton et al., 2025: [Astrocyte and oligodendrocyte pathology in Alzheimer's disease](https://pmc.ncbi.nlm.nih.gov/articles/PMC12047399/)

Publication value:

- **Exploratory and worth mentioning**
- **Not core**

### 8. Reactive hotspots as a novel microenvironmental hypothesis

One of the most interesting outputs from the runs was the repeated emergence of what can reasonably be described as **reactive hotspots** or **reactive focal microenvironments**.

This was most explicit in `run_20260420_215745_90b105`, which produced late accepted features such as:

- `ca1_reactive_hotspot_pyramidal_depletion`
- `ca1_reactive_hotspot_oligodendrocyte_depletion`
- `ca1_reactive_hotspot_corpora_amylacea_enrichment__reactive_minus_astro_hotspot_local_ca_fraction_r70px_k3`
- `ca1_reactive_hotspot_lymphocyte_enrichment__reactive_minus_astro_hotspot_local_lymphocyte_fraction_r70px_k3`

Even outside that run, the same structural idea appeared in other runs under slightly different names:

- local peripyramidal reactive cuffs
- reactive crowding
- reactive local purity
- reactive neuron exclusion
- reactive immune-admixed cuffs
- reactive astrocyte + lymphocyte triad/contact features

Interpretation:

- the system was not only discovering that reactive astrocytes matter
- it was increasingly localizing the signal to **reactive focal niches** and then asking what other cells or debris compartments co-occur there

That is scientifically interesting because it shifts the interpretation from:

- diffuse glial burden

to:

- **microenvironmental organization of pathology**

In other words, the runs suggest that the predictive signal may reside not simply in how many reactive astrocytes exist, but in whether they form **densified hotspot-like niches** characterized by:

- local pyramidal depletion or isolation
- immune/lymphocyte admixture
- corpora amylacea accumulation
- sometimes oligodendrocyte depletion

Novelty judgment:

- The component biology is not novel in isolation.
- The **hotspot-centered quantitative formulation** is more unusual and is one of the clearest novelty candidates produced by these runs.

Caveat:

- The clearest hotspot-heavy run, `run_20260420_215745_90b105`, contains late accepted members with `coverage_ratio = 0.8`, so this should be framed as **novel and hypothesis-generating**, not fully settled.

Publication value:

- **Yes, worth presenting as one of the most interesting novel extensions**
- **Best described as an exploratory microenvironmental hypothesis**

## What Is Publication-Worthy

### Strongest publication-worthy statement

Across five 10-round autoresearch runs, the most reproducible signal was a **CA1-centered reactive astrocyte niche phenotype** in which:

- reactive astrocyte burden/enrichment,
- local peripyramidal reactive organization,
- and pyramidal vulnerability / exclusion

jointly formed the strongest predictive panel for memory decline.

This is publication-worthy because:

- it replicated across independent runs
- it generalized at the **panel level** on held-out test
- the aggregate train/test curve remained coherent across independent runs
- the finding is biologically interpretable

### Strongest novel-ish claim that is still defensible

The system repeatedly converged on **peripyramidal reactive-astrocyte cuffs / crowded niches / neighborhood fractions** rather than only on gross astrocyte burden.

That is the clearest place where the runs appear to add something more specific than standard AD pathology descriptions.

Closely related to that, the most interesting **microenvironmental** extension is that some runs further sharpened the signal into **reactive hotspots** with local pyramidal depletion, immune enrichment, corpora amylacea enrichment, and sometimes oligodendrocyte depletion.

### Strongest exploratory extension

Immune / lymphocyte-associated reactive niches appeared repeatedly as later panel refinements. These are scientifically interesting, but because their standalone significance is mixed, they should be presented as **complementary panel features** rather than uniformly validated biomarkers.

## What Should Not Be Overstated

- The runs did **not** discover an entirely new AD mechanism.
- The runs did **not** prove that every late panel member is a strong standalone biomarker.
- The strongest claims should remain at the **panel** and **convergent motif** level, not at the level of every single round.
- `run_20260420_215745_90b105` should be interpreted with caution for its late `coverage_ratio = 0.8` members.

## Practical Publication Framing

If this material were written up, the most defensible framing would be:

1. **Confirmatory result**
   - CA1 pyramidal vulnerability and reactive astrogliosis are the dominant reproducible histopathologic signal associated with memory decline.

2. **More specific panel result**
   - The predictive signal is sharpened when these features are recast as **local CA1 reactive-astrocyte niches** around vulnerable pyramidal neurons.

3. **Exploratory extension**
   - Some of the most promising later refinements involve **reactive hotspot ecology**, especially immune / lymphocyte admixture, with corpora amylacea and oligodendrocyte gap features as secondary hypotheses worth follow-up.

## One-Sentence Summary

Across the five autoresearch runs, the most publication-worthy finding was not merely that reactive astrocytes matter in AD, but that **CA1 peripyramidal reactive-astrocyte microenvironments repeatedly emerged as the most predictive and biologically coherent panel motif, with hotspot-like immune-, corpora-, and oligodendroglia-associated niche refinements representing the most interesting novel exploratory extension.**
