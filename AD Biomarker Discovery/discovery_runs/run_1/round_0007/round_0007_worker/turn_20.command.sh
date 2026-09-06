python - <<'PY'
from pathlib import Path
path = Path('/scratch/result.py')
text = path.read_text()
old = '''    if worst_ids:
        worst_meta = analyzable.loc[analyzable["donor_id"].isin(worst_ids), ["donor_id", "cognitive_status", "braak_numeric", "cerad_ordinal", best_col]]
        shared_cog = ", ".join(sorted(worst_meta["cognitive_status"].dropna().astype(str).unique().tolist()))
        shared_braak = ", ".join(sorted({str(int(x)) for x in worst_meta["braak_numeric"].dropna().tolist()}))
        shared_cerad = ", ".join(sorted({str(int(x)) for x in worst_meta["cerad_ordinal"].dropna().tolist()}))
        error_pattern = (
            f"The largest raw LOO errors were {', '.join(worst_ids)}; these donors span cognitive states [{shared_cog}] "
            f"with Braak [{shared_braak}] and CERAD [{shared_cerad}], suggesting the niche signal misses some pathology severity heterogeneity."
        )
'''
new = '''    if worst_ids:
        meta_cols = [col for col in ["donor_id", "cognitive_status", "braak_numeric", "cerad_ordinal", best_col] if col in analyzable.columns]
        worst_meta = analyzable.loc[analyzable["donor_id"].isin(worst_ids), meta_cols]
        shared_cog = ", ".join(sorted(worst_meta["cognitive_status"].dropna().astype(str).unique().tolist())) if "cognitive_status" in worst_meta.columns else "unknown"
        shared_braak = ", ".join(sorted({str(int(x)) for x in worst_meta["braak_numeric"].dropna().tolist()})) if "braak_numeric" in worst_meta.columns else "unknown"
        shared_cerad = ", ".join(sorted({str(int(x)) for x in worst_meta["cerad_ordinal"].dropna().tolist()})) if "cerad_ordinal" in worst_meta.columns else "unknown"
        error_pattern = (
            f"The largest raw LOO errors were {', '.join(worst_ids)}; these donors span cognitive states [{shared_cog}] "
            f"with Braak [{shared_braak}] and CERAD [{shared_cerad}], suggesting the niche signal misses some pathology severity heterogeneity."
        )
'''
if old not in text:
    raise SystemExit('target block not found')
path.write_text(text.replace(old, new))
PY
python -m py_compile /scratch/result.py
