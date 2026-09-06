python - <<'PY'
from pathlib import Path
p=Path('/scratch/result.py')
text=p.read_text()
text=text.replace(
'    stronger = "higher" if best["partial_r"] < 0 else "lower"\n',
'    association = "worse" if best["partial_r"] < 0 else "less severe"\n'
)
old = """   - The best signal came from a neuron-centered CA1 niche: the fraction of CA1 pyramidal neurons that simultaneously sit within 30 um of a CA1 reactive astrocyte and within {int(best['lymphocyte_radius_um'])} um of a CA1 lymphocyte. Because worse memory decline corresponds to more negative `slope_zmem0`, the negative partial r means donors with more of this inflammatory perineuronal triad tended to have worse memory decline after confound adjustment.\n"""
new = """   - The best signal came from a neuron-centered CA1 niche: the fraction of CA1 pyramidal neurons that simultaneously sit within 30 um of a CA1 reactive astrocyte and within {int(best['lymphocyte_radius_um'])} um of a CA1 lymphocyte. Because worse memory decline corresponds to more negative `slope_zmem0`, the positive partial r here means donors with more of this perineuronal triad tended to show less severe memory decline after confound adjustment.\n"""
text=text.replace(old,new)
old = """   - The nearby radius alternatives were weaker. The tighter 15 um version appears too sparse, likely missing true immune-adjacent niches when lymphocytes are rare in CA1. The broader 25 um version dilutes locality by counting looser immune proximity that is less specifically perineuronal.\n"""
new = """   - The nearby radius alternatives were weaker. The tighter 15 um version appears too sparse, likely missing true immune-adjacent niches when lymphocytes are rare in CA1. The intermediate 20 um version retained some signal but was still weaker than 25 um, suggesting that this niche only becomes detectable when the lymphocyte neighborhood is allowed to be slightly broader.\n"""
text=text.replace(old,new)
old = """The winning approach is biologically coherent because it sharpens prior astrocyte-only CA1 summaries into a more specific inflammatory niche anchored on the neuron itself. It beat the nearby alternatives because it balances rarity and specificity: 15 um is too strict for a sparse lymphocyte population, while 25 um relaxes the niche enough to mix true local immune-gliosis contact with background immune presence. This family is close to the current panel's CA1 astrocyte and lymphocyte-contact members, so it may add only modestly new information, but it is more interpretable as a neuron-exposure burden.\n"""
new = """The winning approach is biologically coherent because it sharpens prior astrocyte-only CA1 summaries into a more specific inflammatory niche anchored on the neuron itself. It beat the nearby alternatives because 15 um is too strict for a sparse lymphocyte population, whereas 25 um appears to recover a still-local but less brittle version of the same niche. This family is close to the current panel's CA1 astrocyte and lymphocyte-contact members, so it may add only modestly new information, but it is more interpretable as a neuron-exposure burden.\n"""
text=text.replace(old,new)
old = """The signal seems to mean biologically that donors with worse memory decline have a {stronger} burden of CA1 pyramidal neurons embedded in a local inflammatory gliosis niche.  \n"""
new = """The signal seems to mean biologically that donors with {association} memory decline have a greater burden of CA1 pyramidal neurons embedded in a local reactive-astrocyte plus lymphocyte niche.  \n"""
text=text.replace(old,new)
p.write_text(text)
PY
python -m py_compile /scratch/result.py
