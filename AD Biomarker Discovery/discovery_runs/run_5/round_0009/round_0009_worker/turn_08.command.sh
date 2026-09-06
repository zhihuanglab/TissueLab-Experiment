python - <<'PY'
from pathlib import Path
path = Path('/scratch/result.py')
text = path.read_text()
text = text.replace('import json\nimport math\nimport sys\n', 'import json\nimport math\nimport sys\nimport warnings\n')
text = text.replace(
"""from scipy.spatial import cKDTree

# Make shared helpers importable in the sandbox runtime.
""",
"""from scipy.spatial import cKDTree

warnings.filterwarnings("ignore", message=".*not recognized as a component of a Zarr hierarchy.*")

# Make shared helpers importable in the sandbox runtime.
""")
old = """    ranked = sorted(
        variation_results,
        key=lambda x: (
            -999.0 if not math.isfinite(x["selection_score"]) else -x["selection_score"],
            -999.0 if not math.isfinite(x["adjusted_score"]) else -x["adjusted_score"],
            -999.0 if not math.isfinite(x["loo_predictive_r"]) else -x["loo_predictive_r"],
        ),
    )
"""
new = """    def _rank_key(item: dict[str, Any]) -> tuple[float, float, float, float]:
        finite_sel = math.isfinite(item["selection_score"])
        finite_adj = math.isfinite(item["adjusted_score"])
        finite_loo = math.isfinite(item["loo_predictive_r"])
        return (
            1.0 if finite_sel else 0.0,
            item["selection_score"] if finite_sel else -np.inf,
            item["adjusted_score"] if finite_adj else -np.inf,
            item["loo_predictive_r"] if finite_loo else -np.inf,
        )

    ranked = sorted(variation_results, key=_rank_key, reverse=True)
"""
if old not in text:
    raise SystemExit('ranking block not found')
text = text.replace(old, new)
path.write_text(text)
PY

python /scratch/result.py | tee /scratch/result_stdout.txt
