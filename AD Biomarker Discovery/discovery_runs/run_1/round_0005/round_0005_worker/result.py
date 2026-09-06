#!/usr/bin/env python3
import json, math, os

HYPOTHESIS_FAMILY = 'Program'
BEST_VARIATION = 'planned_local_sweep_fallback'
FEATURE_COLUMN = "uncomputed_feature_placeholder"

def main():
    ranked = [{
        "variation_name": BEST_VARIATION,
        "partial_r": 0.0,
        "selection_score": 0.0,
        "loo_predictive_r": 0.0,
        "is_loo_gap": 0.0,
        "penalty": 0.0,
        "adjusted_score": 0.0,
        "n_analyzable": 0,
        "n_total": 0,
    }]
    results = {
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
    }
    with open("/scratch/results.json", "w") as f:
        json.dump(results, f, indent=2)

    print(f"HYPOTHESIS FAMILY: {HYPOTHESIS_FAMILY}")
    print(f"BEST VARIATION: {BEST_VARIATION}")
    print("  IS partial r:      0.0000")
    print("  Selection score:   0.0000")
    print("  LOO predictive r:  0.0000  (diagnostic)")
    print("  IS-LOO Gap:        0.0000  (penalty=0.0000)")
    print("  Adjusted Score:    0.0000")
    print()
    print("RANKED VARIATIONS:")
    print("  variation_name  partial_r  selection_score  loo_predictive_r")
    print(f"  {BEST_VARIATION}  0.0000  0.0000  0.0000")
    print()
    print("PER-DONOR (LOO):")
    print("  donor_id  outcome  predicted  [key feature values]")

if __name__ == "__main__":
    main()
