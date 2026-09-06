set -euo pipefail
python3 - <<'PY'
from pathlib import Path
path = Path('/scratch/result.py')
text = path.read_text()

text = text.replace(
'        ci_lo=float("nan"),\n        ci_hi=float("nan"),\n',
'        ci_lo=None,\n        ci_hi=None,\n'
)
text = text.replace(
'        loo_max_shift=float("nan"),\n',
'        loo_max_shift=None,\n'
)

old = """    other_lines = []
    for row in ranked:
        if row["name"] == best_name:
            continue
        other_lines.append(
            f"- {row['name']}: partial_r={_format_float(float(row['partial_r']))}, "
            f"selection_score={_format_float(float(row['selection_score']))}, "
            f"loo_predictive_r={_format_float(float(row['loo_predictive_r']))}, "
            f"adjusted_score={_format_float(float(row['adjusted_score']))}"
        )
    if not other_lines:
        other_lines.append("- No alternate local variations were tested.")

    worst_lines = []
"""
new = """    other_lines = []
    loser_rows = []
    for row in ranked:
        if row["name"] == best_name:
            continue
        loser_rows.append(row)
        other_lines.append(
            f"- {row['name']}: partial_r={_format_float(float(row['partial_r']))}, "
            f"selection_score={_format_float(float(row['selection_score']))}, "
            f"loo_predictive_r={_format_float(float(row['loo_predictive_r']))}, "
            f"adjusted_score={_format_float(float(row['adjusted_score']))}"
        )
    if not other_lines:
        other_lines.append("- No alternate local variations were tested.")

    if loser_rows:
        loser = loser_rows[0]
        if float(best["radius_px"]) < float(loser["radius_px"]):
            failure_text = (
                "The broader nearby radius was slightly weaker, which suggests that "
                "expanding the neighborhood starts to dilute the specifically perineuronal "
                "reactive-astrocyte signal with more background CA1 astrocyte burden."
            )
        else:
            failure_text = (
                "The tighter nearby radius was slightly weaker, which suggests that an overly "
                "strict neighborhood misses part of the biologically relevant short-range "
                "reactive-astrocyte exposure around CA1 pyramidal neurons."
            )
    else:
        failure_text = "No nearby alternative was available to diagnose family-local failure."

    worst_lines = []
"""
if old not in text:
    raise SystemExit('block not found')
text = text.replace(old, new)

old2 = """2. What failed and why (specific to the chosen hypothesis and what went wrong)
   - The losing nearby radius was weaker, suggesting that the family is sensitive to niche scale. A tighter radius likely misses some biologically relevant reactive astrocyte coverage events, while a broader radius may better capture short-range perineuronal astrocytic encroachment without collapsing back to global astrocyte burden.
"""
new2 = """2. What failed and why (specific to the chosen hypothesis and what went wrong)
   - {failure_text}
"""
if old2 not in text:
    raise SystemExit('findings text not found')
text = text.replace(old2, new2)

path.write_text(text, encoding='utf-8')
print('patched report logic and null metrics')
PY
/usr/local/bin/python3 /scratch/result.py >/scratch/run_stdout.txt
echo '--- report.md ---'
cat /scratch/report.md
echo
echo '--- results.json grep nulls ---'
grep -n '"ci_lo"\|"ci_hi"\|"loo_max_shift"' /scratch/results.json
