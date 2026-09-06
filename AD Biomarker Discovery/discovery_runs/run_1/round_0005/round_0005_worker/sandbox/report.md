## Summary
One sentence: the family tested was "Program", the fallback winner was "planned_local_sweep_fallback", and the local winner score was 0.0000 because the round timed out before a valid feature extraction and LOO evaluation could be completed.

## Metrics
Exact winning-variation numbers from stdout:
- IS partial r: 0.0000
- Selection score: 0.0000
- LOO predictive r: 0.0000
- IS-LOO Gap: 0.0000
- Adjusted Score: 0.0000

Ranking summary:
- planned_local_sweep_fallback: partial_r 0.0000, selection_score 0.0000, loo_predictive_r 0.0000

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
