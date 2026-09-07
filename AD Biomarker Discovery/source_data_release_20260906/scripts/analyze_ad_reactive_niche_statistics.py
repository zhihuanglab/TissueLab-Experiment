#!/usr/bin/env python3
"""Comprehensive statistics for the rebuilt CA1 reactive-astrocyte niche."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr, t as student_t
from statsmodels.stats.multitest import multipletests


COVARIATES = ["max_age_vis", "sex", "braak_numeric", "cerad_ordinal"]
OUTCOME = "slope_zmem0"
BULK = "ca1_reactive_fraction"

SPATIAL_SPECS = [
    ("reactive_hit_over_astroglial_hit_35um", -1, "Reactive/astroglial-neighbored pyramidal ratio, 35 um", "core"),
    ("peripyramidal_lineage_reactivity_35um", -1, "Peripyramidal lineage reactivity, 35 um", "core"),
    ("peripyramidal_lineage_reactivity_30um", -1, "Peripyramidal lineage reactivity, 30 um", "core"),
    ("pyr_ge1_reactive_40um", -1, ">=1 reactive astrocyte within 40 um", "core"),
    ("pyr_ge1_reactive_35um", -1, ">=1 reactive astrocyte within 35 um", "core"),
    ("pyr_ge1_reactive_30um", -1, ">=1 reactive astrocyte within 30 um", "core"),
    ("pyr_supermajority_reactive_cuff_35um", -1, "Supermajority-reactive cuff, 35 um", "core"),
    ("pyr_ge2_reactive_35um", -1, ">=2 reactive astrocytes within 35 um", "core"),
    ("pyr_ge3_reactive_30um", -1, ">=3 reactive astrocytes within 30 um", "core"),
    ("reactive_local_purity_60um", -1, "Reactive local purity, 60 um", "refinement"),
    ("reactive_pyramidal_desert_t20_35um", -1, "Reactive-centered pyramidal desert", "refinement"),
    ("pyr_ge1_reactive_literal_6px", 1, "Literal 6-px proximity (non-transportable)", "degenerate"),
]

DUPLICATE_FORMULA = "pyr_reactive_isolated_k5_30um"
DUPLICATE_OF = "pyr_ge1_reactive_30um"


def sex_numeric(series: pd.Series) -> np.ndarray:
    return (
        series.astype(str)
        .str.strip()
        .str.lower()
        .map({"female": 0.0, "f": 0.0, "male": 1.0, "m": 1.0})
        .to_numpy(dtype=float)
    )


def covariate_matrix(frame: pd.DataFrame) -> np.ndarray:
    return np.column_stack(
        [
            np.ones(len(frame)),
            frame["max_age_vis"].to_numpy(dtype=float),
            sex_numeric(frame["sex"]),
            frame["braak_numeric"].to_numpy(dtype=float),
            frame["cerad_ordinal"].to_numpy(dtype=float),
        ]
    )


def residualize(values: np.ndarray, covariates: np.ndarray) -> np.ndarray:
    return values - covariates @ np.linalg.lstsq(covariates, values, rcond=None)[0]


def partial_r(feature: np.ndarray, outcome: np.ndarray, covariates: np.ndarray) -> float:
    x = residualize(feature, covariates)
    y = residualize(outcome, covariates)
    if np.std(x) <= 1e-15 or np.std(y) <= 1e-15:
        return float("nan")
    return float(np.corrcoef(x, y)[0, 1])


def residual_spearman(feature: np.ndarray, outcome: np.ndarray, covariates: np.ndarray) -> float:
    x = residualize(feature, covariates)
    y = residualize(outcome, covariates)
    return float(spearmanr(x, y).statistic)


def bootstrap_metrics(
    feature: np.ndarray,
    outcome: np.ndarray,
    covariates: np.ndarray,
    *,
    expected_sign: int,
    n_bootstrap: int,
    seed: int,
) -> dict[str, float]:
    observed = partial_r(feature, outcome, covariates)
    if not np.isfinite(observed):
        return {
            "bootstrap_ci_low": np.nan,
            "bootstrap_ci_high": np.nan,
            "bootstrap_median": np.nan,
            "bootstrap_declared_direction_fraction": np.nan,
        }
    rng = np.random.default_rng(seed)
    values = []
    for _ in range(n_bootstrap):
        index = rng.integers(0, len(feature), len(feature))
        value = partial_r(feature[index], outcome[index], covariates[index])
        if np.isfinite(value):
            values.append(value)
    values = np.asarray(values)
    low, high = np.quantile(values, [0.025, 0.975])
    return {
        "bootstrap_ci_low": float(low),
        "bootstrap_ci_high": float(high),
        "bootstrap_median": float(np.median(values)),
        "bootstrap_declared_direction_fraction": float(np.mean(expected_sign * values > 0)),
    }


def lodo_metrics(
    frame: pd.DataFrame,
    feature_name: str,
    expected_sign: int,
) -> dict[str, object]:
    full_covariates = covariate_matrix(frame)
    full_r = partial_r(
        frame[feature_name].to_numpy(float),
        frame[OUTCOME].to_numpy(float),
        full_covariates,
    )
    values = []
    for index, row in frame.reset_index(drop=True).iterrows():
        reduced = frame.drop(frame.index[index]).reset_index(drop=True)
        value = partial_r(
            reduced[feature_name].to_numpy(float),
            reduced[OUTCOME].to_numpy(float),
            covariate_matrix(reduced),
        )
        values.append((str(row["donor_id"]), value))
    finite = [(donor, value) for donor, value in values if np.isfinite(value)]
    most_influential = max(finite, key=lambda item: abs(item[1] - full_r))
    correlations = np.asarray([value for _, value in finite])
    return {
        "formula": feature_name,
        "n": len(frame),
        "full_partial_r": full_r,
        "lodo_min_r": float(np.min(correlations)),
        "lodo_max_r": float(np.max(correlations)),
        "lodo_median_r": float(np.median(correlations)),
        "lodo_declared_direction_fraction": float(np.mean(expected_sign * correlations > 0)),
        "most_influential_donor": most_influential[0],
        "max_abs_lodo_change": float(abs(most_influential[1] - full_r)),
    }


def freedman_lane_matrix(
    frame: pd.DataFrame,
    feature_names: list[str],
    expected_signs: np.ndarray,
    *,
    n_permutations: int,
    seed: int,
) -> tuple[pd.DataFrame, dict[str, float]]:
    covariates = covariate_matrix(frame)
    projection = np.eye(len(frame)) - covariates @ np.linalg.pinv(covariates)
    outcome = frame[OUTCOME].to_numpy(float)
    outcome_fit = covariates @ np.linalg.lstsq(covariates, outcome, rcond=None)[0]
    outcome_residual = outcome - outcome_fit
    features = projection @ frame[feature_names].to_numpy(float)
    norms = np.linalg.norm(features, axis=0)
    if np.any(norms <= 1e-15):
        raise RuntimeError("A feature in the permutation family is constant after residualization")
    standardized_features = features / norms
    observed_outcome = projection @ outcome
    observed = observed_outcome / np.linalg.norm(observed_outcome) @ standardized_features
    observed_aligned = observed * expected_signs

    rng = np.random.default_rng(seed)
    raw_two_count = np.zeros(len(feature_names), dtype=np.int64)
    raw_directional_count = np.zeros(len(feature_names), dtype=np.int64)
    max_abs_null = np.empty(n_permutations)
    max_aligned_null = np.empty(n_permutations)
    batch_size = 1000
    for start in range(0, n_permutations, batch_size):
        batch = min(batch_size, n_permutations - start)
        permutations = np.asarray([rng.permutation(len(frame)) for _ in range(batch)])
        permuted_outcome = outcome_fit[None, :] + outcome_residual[permutations]
        residuals = permuted_outcome @ projection.T
        residuals /= np.linalg.norm(residuals, axis=1, keepdims=True)
        correlations = residuals @ standardized_features
        absolute = np.abs(correlations)
        aligned = correlations * expected_signs[None, :]
        raw_two_count += np.sum(absolute >= np.abs(observed)[None, :], axis=0)
        raw_directional_count += np.sum(aligned >= observed_aligned[None, :], axis=0)
        max_abs_null[start : start + batch] = np.max(absolute, axis=1)
        max_aligned_null[start : start + batch] = np.max(aligned, axis=1)

    denominator = n_permutations + 1
    raw_two = (raw_two_count + 1) / denominator
    raw_directional = (raw_directional_count + 1) / denominator
    max_two = np.asarray(
        [(1 + np.sum(max_abs_null >= abs(value))) / denominator for value in observed]
    )
    max_directional = np.asarray(
        [(1 + np.sum(max_aligned_null >= value)) / denominator for value in observed_aligned]
    )
    result = pd.DataFrame(
        {
            "formula": feature_names,
            "partial_r": observed,
            "declared_direction": expected_signs,
            "direction_aligned_r": observed_aligned,
            "permutation_p_two_sided": raw_two,
            "permutation_p_directional": raw_directional,
            "maxT_fwer_p_two_sided": max_two,
            "maxT_fwer_p_directional": max_directional,
        }
    )
    result["bh_fdr_p_two_sided"] = multipletests(raw_two, method="fdr_bh")[1]
    result["bh_fdr_p_directional"] = multipletests(raw_directional, method="fdr_bh")[1]
    omnibus = {
        "n": len(frame),
        "n_formulas": len(feature_names),
        "observed_max_abs_r": float(np.max(np.abs(observed))),
        "maxT_omnibus_p_two_sided": float(
            (1 + np.sum(max_abs_null >= np.max(np.abs(observed)))) / denominator
        ),
        "observed_max_direction_aligned_r": float(np.max(observed_aligned)),
        "maxT_omnibus_p_directional": float(
            (1 + np.sum(max_aligned_null >= np.max(observed_aligned))) / denominator
        ),
    }
    return result, omnibus


def build_current_frames(args: argparse.Namespace) -> tuple[pd.DataFrame, pd.DataFrame]:
    train = pd.read_csv(args.spatial_train)
    test = pd.read_csv(args.spatial_test)
    bulk_train = pd.read_csv(args.bulk_train)[["donor_id", BULK]]
    bulk_test = pd.read_csv(args.bulk_test)[["donor_id", BULK]]
    train = train.merge(bulk_train, on="donor_id", validate="one_to_one")
    test = test.merge(bulk_test, on="donor_id", validate="one_to_one")
    return train, test


def summarize_formula(
    frame: pd.DataFrame,
    formula: str,
    expected_sign: int,
    *,
    split: str,
    n_bootstrap: int,
    seed: int,
) -> dict[str, object]:
    columns = ["donor_id", formula, OUTCOME, *COVARIATES]
    subset = frame[columns].replace([np.inf, -np.inf], np.nan).dropna().reset_index(drop=True)
    feature = subset[formula].to_numpy(float)
    outcome = subset[OUTCOME].to_numpy(float)
    covariates = covariate_matrix(subset)
    record: dict[str, object] = {
        "split": split,
        "formula": formula,
        "n": len(subset),
        "partial_r": partial_r(feature, outcome, covariates),
        "residual_spearman_r": residual_spearman(feature, outcome, covariates),
        "declared_direction": expected_sign,
    }
    record.update(
        bootstrap_metrics(
            feature,
            outcome,
            covariates,
            expected_sign=expected_sign,
            n_bootstrap=n_bootstrap,
            seed=seed,
        )
    )
    return record


def historical_full_search(args: argparse.Namespace, output_dir: Path) -> dict[str, object]:
    historical = pd.read_csv(args.historical_rounds)
    evaluable = historical.loc[historical["test_partial_r"].notna()].copy()
    unique = evaluable.drop_duplicates(["train_partial_r", "test_partial_r"]).copy()
    for frame, name in [(evaluable, "historical_48_evaluable"), (unique, "historical_42_unique_effect_profiles")]:
        degrees_freedom = 34 - len(COVARIATES) - 2
        absolute = frame["test_partial_r"].abs().to_numpy(float)
        statistic = absolute * np.sqrt(degrees_freedom / (1 - absolute**2))
        frame["analytic_p_two_sided"] = 2 * student_t.sf(statistic, degrees_freedom)
        frame["bh_fdr_p_two_sided"] = multipletests(
            frame["analytic_p_two_sided"], method="fdr_bh"
        )[1]
        frame["holm_p_two_sided"] = multipletests(
            frame["analytic_p_two_sided"], method="holm"
        )[1]
        frame.to_csv(output_dir / f"{name}.csv", index=False)
    significant = unique.loc[unique["bh_fdr_p_two_sided"] < 0.05]
    return {
        "round_representatives": int(len(historical)),
        "evaluable_test_associations": int(len(evaluable)),
        "unique_exact_effect_profiles": int(len(unique)),
        "bh_significant_unique_profiles": int(len(significant)),
        "holm_significant_unique_profiles": int((unique["holm_p_two_sided"] < 0.05).sum()),
        "significant_families": significant["family"].value_counts().to_dict(),
        "caveat": "Analytic partial-correlation p-values; exact-profile deduplication is a sensitivity audit, not the primary permutation inference.",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--spatial-train", type=Path, required=True)
    parser.add_argument("--spatial-test", type=Path, required=True)
    parser.add_argument("--bulk-train", type=Path, required=True)
    parser.add_argument("--bulk-test", type=Path, required=True)
    parser.add_argument("--historical-rounds", type=Path, required=True)
    parser.add_argument("--april-transfer", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--n-permutations", type=int, default=100_000)
    parser.add_argument("--n-bootstrap", type=int, default=10_000)
    parser.add_argument("--seed", type=int, default=20260825)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    train, test = build_current_frames(args)
    spatial_names = [spec[0] for spec in SPATIAL_SPECS]
    spatial_signs = np.asarray([spec[1] for spec in SPATIAL_SPECS], dtype=int)
    labels = {spec[0]: spec[2] for spec in SPATIAL_SPECS}
    roles = {spec[0]: spec[3] for spec in SPATIAL_SPECS}

    duplicate_train_diff = float(
        np.max(np.abs(train[DUPLICATE_FORMULA] - train[DUPLICATE_OF]))
    )
    duplicate_test_diff = float(
        np.max(np.abs(test[DUPLICATE_FORMULA] - test[DUPLICATE_OF]))
    )
    pd.DataFrame(
        [
            {
                "formula": DUPLICATE_FORMULA,
                "duplicate_of": DUPLICATE_OF,
                "max_abs_difference_train": duplicate_train_diff,
                "max_abs_difference_test": duplicate_test_diff,
                "decision": "Retain in formula inventory but collapse for distinct-measurement multiplicity.",
            }
        ]
    ).to_csv(args.output_dir / "duplicate_formula_audit.csv", index=False)

    formula_rows = []
    all_specs = [(BULK, -1, "Bulk CA1 reactive-astrocyte fraction", "bulk"), *SPATIAL_SPECS]
    for split, frame in [("Train", train), ("Test", test)]:
        for offset, (formula, expected_sign, label, role) in enumerate(all_specs):
            record = summarize_formula(
                frame,
                formula,
                expected_sign,
                split=split,
                n_bootstrap=args.n_bootstrap,
                seed=args.seed + (0 if split == "Train" else 10_000) + offset,
            )
            record["label"] = label
            record["role"] = role
            formula_rows.append(record)
    formula_statistics = pd.DataFrame(formula_rows)

    permutation_table, spatial_omnibus = freedman_lane_matrix(
        test,
        spatial_names,
        spatial_signs,
        n_permutations=args.n_permutations,
        seed=args.seed,
    )
    permutation_table["label"] = permutation_table["formula"].map(labels)
    permutation_table["role"] = permutation_table["formula"].map(roles)
    test_statistics = formula_statistics.loc[formula_statistics["split"].eq("Test")]
    permutation_table = permutation_table.merge(
        test_statistics.drop(columns=["split", "label", "role", "partial_r", "declared_direction"]),
        on="formula",
        validate="one_to_one",
    )
    permutation_table.to_csv(args.output_dir / "spatial_formula_statistics.csv", index=False)
    formula_statistics.to_csv(args.output_dir / "train_test_formula_statistics.csv", index=False)

    bulk_perm, bulk_omnibus = freedman_lane_matrix(
        test,
        [BULK],
        np.asarray([-1]),
        n_permutations=args.n_permutations,
        seed=args.seed + 1_000_000,
    )
    component_p = np.asarray(
        [
            float(bulk_perm.iloc[0]["permutation_p_two_sided"]),
            float(spatial_omnibus["maxT_omnibus_p_two_sided"]),
        ]
    )
    component_directional_p = np.asarray(
        [
            float(bulk_perm.iloc[0]["permutation_p_directional"]),
            float(spatial_omnibus["maxT_omnibus_p_directional"]),
        ]
    )
    component_table = pd.DataFrame(
        [
            {
                "component": "Bulk CA1 reactive-astrocyte enrichment",
                "n": len(test),
                "effect_statistic": float(bulk_perm.iloc[0]["partial_r"]),
                "effect_statistic_definition": "partial Pearson r",
                "permutation_p_two_sided": component_p[0],
                "permutation_p_directional": component_directional_p[0],
            },
            {
                "component": "Peripyramidal reactive-astrocyte spatial family",
                "n": len(test),
                "effect_statistic": spatial_omnibus["observed_max_abs_r"],
                "effect_statistic_definition": "maximum absolute partial r across 12 distinct measurements",
                "permutation_p_two_sided": component_p[1],
                "permutation_p_directional": component_directional_p[1],
            },
        ]
    )
    component_table["holm_p_two_sided_across_two_components"] = multipletests(
        component_p, method="holm"
    )[1]
    component_table["holm_p_directional_across_two_components"] = multipletests(
        component_directional_p, method="holm"
    )[1]
    component_table.to_csv(args.output_dir / "primary_component_statistics.csv", index=False)

    combined_family_table, combined_family_omnibus = freedman_lane_matrix(
        test,
        [BULK, *spatial_names],
        np.asarray([-1, *spatial_signs]),
        n_permutations=args.n_permutations,
        seed=args.seed + 2_000_000,
    )
    combined_family_table.to_csv(
        args.output_dir / "combined_family_formula_statistics.csv", index=False
    )

    lodo_rows = []
    for formula, expected_sign, _, _ in all_specs:
        subset = test[["donor_id", formula, OUTCOME, *COVARIATES]].dropna().reset_index(drop=True)
        lodo_rows.append(lodo_metrics(subset, formula, expected_sign))
    pd.DataFrame(lodo_rows).to_csv(args.output_dir / "test_lodo_stability.csv", index=False)

    # Residualized feature redundancy on held-out donors.
    feature_names = [BULK, *spatial_names]
    residual_features = residualize(
        test[feature_names].to_numpy(float), covariate_matrix(test)
    )
    correlation = pd.DataFrame(
        np.corrcoef(residual_features, rowvar=False),
        index=feature_names,
        columns=feature_names,
    )
    correlation.to_csv(args.output_dir / "test_residualized_feature_correlations.csv")

    transfer = pd.read_csv(args.april_transfer)
    transfer_valid = transfer.replace([np.inf, -np.inf], np.nan).dropna(
        subset=["original_train_r", "original_test_r", "rebuilt_train_r", "rebuilt_test_r"]
    )
    transfer_summary = {
        "n_april_entries": int(len(transfer)),
        "n_entries_with_complete_effects": int(len(transfer_valid)),
        "original_train_vs_rebuilt_test_pearson": float(
            pearsonr(transfer_valid["original_train_r"], transfer_valid["rebuilt_test_r"]).statistic
        ),
        "original_test_vs_rebuilt_test_pearson": float(
            pearsonr(transfer_valid["original_test_r"], transfer_valid["rebuilt_test_r"]).statistic
        ),
        "original_test_vs_rebuilt_test_spearman": float(
            spearmanr(transfer_valid["original_test_r"], transfer_valid["rebuilt_test_r"]).statistic
        ),
        "rebuilt_test_direction_preserved_count": int(transfer["test_direction_preserved"].sum()),
    }

    historical_summary = historical_full_search(args, args.output_dir)
    sample_sizes = pd.DataFrame(
        [
            ["Discovery associations", 35, "Hypothesis generation; selection-inflated"],
            ["Held-out component and family inference", 34, "Primary inferential unit"],
            ["Five agent runs", 35, "Search reproducibility on the same donors; not n=175"],
            ["Rebuilt-stack replay", 34, "Measurement robustness on the same Test donors"],
            ["Pooled Train + Test", 69, "Descriptive only for outcome association"],
        ],
        columns=["analysis", "n_donors", "role"],
    )
    sample_sizes.to_csv(args.output_dir / "sample_size_roles.csv", index=False)

    summary = {
        "generated_from": {
            "spatial_train": str(args.spatial_train),
            "spatial_test": str(args.spatial_test),
            "bulk_train": str(args.bulk_train),
            "bulk_test": str(args.bulk_test),
            "historical_rounds": str(args.historical_rounds),
        },
        "settings": {
            "covariates": COVARIATES,
            "outcome": OUTCOME,
            "n_permutations": args.n_permutations,
            "n_bootstrap": args.n_bootstrap,
            "seed": args.seed,
            "primary_sidedness": "two-sided",
        },
        "spatial_omnibus": spatial_omnibus,
        "bulk_omnibus": bulk_omnibus,
        "combined_reactive_astrocyte_family_omnibus": combined_family_omnibus,
        "duplicate_formula": {
            "formula": DUPLICATE_FORMULA,
            "duplicate_of": DUPLICATE_OF,
            "max_abs_difference_train": duplicate_train_diff,
            "max_abs_difference_test": duplicate_test_diff,
        },
        "transfer_summary": transfer_summary,
        "historical_full_search_summary": historical_summary,
    }
    (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(component_table.to_string(index=False))
    print("\nSpatial omnibus:", json.dumps(spatial_omnibus, indent=2))
    print("\nHistorical audit:", json.dumps(historical_summary, indent=2))


if __name__ == "__main__":
    main()
