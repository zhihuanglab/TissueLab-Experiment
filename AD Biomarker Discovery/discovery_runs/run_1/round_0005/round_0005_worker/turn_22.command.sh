set -e
python - <<'PY'
import json, os, re, textwrap, subprocess, sys

context_path = "/scratch/context_bundle.json"
program_path = "/scratch/program.md"

family = "unresolved_candidate_family"
variation = "fallback_null_variation"

ctx = {}
if os.path.exists(context_path):
    try:
        with open(context_path) as f:
            ctx = json.load(f)
    except Exception:
        ctx = {}
prog = ""
if os.path.exists(program_path):
    try:
        prog = open(program_path).read()
    except Exception:
        prog = ""

def pick_family():
    # Try a few likely fields first
    for k in ["hypothesis_family", "family", "candidate_family", "round_name", "title", "hypothesis"]:
        v = ctx.get(k)
        if isinstance(v, str) and v.strip():
            return v.strip()[:120]
    # Search nested dict/list strings
    txt = json.dumps(ctx)[:8000]
    m = re.search(r'"(?:hypothesis_family|candidate_family|family|round_name|title)"\s*:\s*"([^"]+)"', txt)
    if m:
        return m.group(1)[:120]
    # Search markdown
    for pat in [
        r'(?im)^#\s+(.+)$',
        r'(?im)^##\s+(.+)$',
        r'(?i)family[:\s]+([^\n]+)',
        r'(?i)candidate[:\s]+([^\n]+)',
        r'(?i)hypothesis[:\s]+([^\n]+)',
    ]:
        m = re.search(pat, prog)
        if m:
            s = m.group(1).strip()
            if s:
                return s[:120]
    return family

family = pick_family()
# derive a slightly more informative variation label if possible
m = re.search(r'(?i)(variation|sweep|compare|test)(.+)', prog)
if m:
    variation = "planned_local_sweep_fallback"
    
result_py = f'''#!/usr/bin/env python3
import json, math, os

HYPOTHESIS_FAMILY = {family!r}
BEST_VARIATION = {variation!r}
FEATURE_COLUMN = "uncomputed_feature_placeholder"

def main():
    ranked = [{{
        "variation_name": BEST_VARIATION,
        "partial_r": 0.0,
        "selection_score": 0.0,
        "loo_predictive_r": 0.0,
        "is_loo_gap": 0.0,
        "penalty": 0.0,
        "adjusted_score": 0.0,
        "n_analyzable": 0,
        "n_total": 0,
    }}]
    results = {{
        "best_variation": BEST_VARIATION,
        "feature_column": FEATURE_COLUMN,
        "partial_r": 0.0,
        "selection_score": 0.0,
        "loo_predictive_r": 0.0,
        "is_loo_gap": 0.0,
        "penalty": 0.0,
        "adjusted_score": 0.0,
        "ranked_variations": ranked,
        "note": "Emergency fallback result written because the round could not be completed before timeout. This preserves replay structure and a canonical target for later replacement."
    }}
    with open("/scratch/results.json", "w") as f:
        json.dump(results, f, indent=2)

    print(f"HYPOTHESIS FAMILY: {{HYPOTHESIS_FAMILY}}")
    print(f"BEST VARIATION: {{BEST_VARIATION}}")
    print("  IS partial r:      0.0000")
    print("  Selection score:   0.0000")
    print("  LOO predictive r:  0.0000  (diagnostic)")
    print("  IS-LOO Gap:        0.0000  (penalty=0.0000)")
    print("  Adjusted Score:    0.0000")
    print()
    print("RANKED VARIATIONS:")
    print("  variation_name  partial_r  selection_score  loo_predictive_r")
    print(f"  {{BEST_VARIATION}}  0.0000  0.0000  0.0000")
    print()
    print("PER-DONOR (LOO):")
    print("  donor_id  outcome  predicted  [key feature values]")

if __name__ == "__main__":
    main()
'''
open("/scratch/result.py","w").write(result_py)
os.chmod("/scratch/result.py", 0o755)

out = subprocess.check_output([sys.executable, "/scratch/result.py"], text=True)

report = f"""## Summary
One sentence: the family tested was "{family}", the fallback winner was "{variation}", and the local winner score was 0.0000 because the round timed out before a valid feature extraction and LOO evaluation could be completed.

## Metrics
Exact winning-variation numbers from stdout:
- IS partial r: 0.0000
- Selection score: 0.0000
- LOO predictive r: 0.0000
- IS-LOO Gap: 0.0000
- Adjusted Score: 0.0000

Ranking summary:
- {variation}: partial_r 0.0000, selection_score 0.0000, loo_predictive_r 0.0000

## Findings
1. What worked and why: the round at least preserves the intended replay contract by writing a canonical `result.py`, `results.json`, and a named best variation so the candidate can be revisited cleanly.
2. What failed and why: no valid donor-level scalar was computed before timeout, so the planned biomarker family could not be evaluated and all metrics are placeholder null values.
3. Error pattern: no donor-level prediction residuals were produced, so no specific donor error pattern can be claimed from this fallback output.

## Rationale
This file is an honest null placeholder rather than an analyzed biomarker. It does not support a biological claim, and it should not be interpreted as evidence for or against the family. Its only purpose is to preserve a deterministic replay target and indicate that the round needs rerun.

## Interpretation
No biological interpretation is justified from this fallback result. The population, niche, feature summary, and observable tissue pattern remain unresolved for this round.

## Next
Re-run this exact round with enough time to finish the planned local sweep, extract the intended donor-level feature once outside CV, and then evaluate the nearby family variations under full LOO as originally specified.
"""
open("/scratch/report.md","w").write(report)
print(out)
PY
ls -l /scratch/result.py /scratch/report.md /scratch/results.json
